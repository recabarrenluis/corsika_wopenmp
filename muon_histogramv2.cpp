#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TClass.h>
#include <TParameter.h>
#include <TString.h>

#include <vector>
#include <iostream>
#include <cmath>
#include <memory>
#include <omp.h>


// change ID particle to extract other particle histogram. E.g. 2,3 for electron/positron
// DAT.root as parsed argument
int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <root_file>\n";
        return 1;
    }

    // runtime computation
    double t0 = omp_get_wtime();

    // these lines are necessary to prevent multithreading crash when reading root files
    std::string filepath = std::string("./data/") + argv[1];
    const char* filename = filepath.c_str();
    ROOT::EnableThreadSafety();

    gROOT->ProcessLine("#include <vector>");

    // warm-up ROOT dictionaries / streamer info
    TFile warm_file(filename, "READ");
    TTree* warm_tree = (TTree*)warm_file.Get("events");

    std::vector<float>* warm_px = nullptr;
    std::vector<float>* warm_py = nullptr;
    std::vector<float>* warm_pz = nullptr;
    std::vector<short>* warm_particle_type = nullptr;

    warm_tree->SetBranchAddress("px", &warm_px);
    warm_tree->SetBranchAddress("py", &warm_py);
    warm_tree->SetBranchAddress("pz", &warm_pz);
    warm_tree->SetBranchAddress("particle_type", &warm_particle_type);
    warm_tree->GetEntry(0);

    warm_file.Close();

    // after initializing ROOT to multithreading well-functioning we open our file
    TFile file(filename, "READ");
    TTree* tree = (TTree*)file.Get("events");

    int nEvents = tree->GetEntries();

    // Read primary info
    auto primary_energy  = (TParameter<float>*)file.Get("primary_energy");
    auto primary_zenith  = (TParameter<float>*)file.Get("primary_zenith");
    auto primary_azimuth = (TParameter<float>*)file.Get("primary_azimuth");
    auto primary_type    = (TParameter<int>*)file.Get("primary_type");

    std::cout << "*** PRIMARY INFO ***\n";
    if (primary_energy)
        std::cout << "Energy  = " << primary_energy->GetVal() << std::endl;
    if (primary_zenith)
        std::cout << "Zenith  = " << primary_zenith->GetVal() << std::endl;
    if (primary_azimuth)
        std::cout << "Azimuth = " << primary_azimuth->GetVal() << std::endl;
    if (primary_type)
        std::cout << "Type    = " << primary_type->GetVal() << std::endl;

    // find energy range (for histogram bin) and total muons
    double global_min = 1e30;
    double global_max = -1e30;
    long long total_muons = 0;

#pragma omp parallel
    {
        TFile f(filename, "READ");
        TTree* t = (TTree*)f.Get("events");

        std::vector<float>* px = nullptr;
        std::vector<float>* py = nullptr;
        std::vector<float>* pz = nullptr;
        std::vector<short>* particle_type = nullptr;

        t->SetBranchAddress("px", &px);
        t->SetBranchAddress("py", &py);
        t->SetBranchAddress("pz", &pz);
        t->SetBranchAddress("particle_type", &particle_type);

        double local_min = 1e30;
        double local_max = -1e30;
        long long local_count = 0;

#pragma omp for schedule(dynamic,1)
        for (int i = 0; i < nEvents; i++)
        {
            t->GetEntry(i);

            int n = px->size();

            for (int j = 0; j < n; j++)
            {
                short pid = particle_type->at(j);

                if (pid == 5 || pid == 6)
                {
                    float e = std::sqrt(
                        px->at(j) * px->at(j) +
                        py->at(j) * py->at(j) +
                        pz->at(j) * pz->at(j)
                    );

                    if (e < local_min) local_min = e;
                    if (e > local_max) local_max = e;

                    local_count++;
                }
            }
        }

#pragma omp critical
        {
            if (local_min < global_min) global_min = local_min;
            if (local_max > global_max) global_max = local_max;
            total_muons += local_count;
        }

        f.Close();
    }

    std::cout << "Total muons = " << total_muons << std::endl;

    int bins = static_cast<int>(sqrt(total_muons));
    std::cout << "Histogram bins = " << bins << std::endl;

    std::vector<double> bin_edges(bins + 1);

    double log_min = log10(global_min);
    double log_max = log10(global_max);

    for (int i = 0; i <= bins; i++)
    {
        bin_edges[i] = pow(10, log_min + (log_max - log_min) * i / bins);
    }

    auto h_muon_energy = std::make_unique<TH1F>(
        "h_muon_energy",
        "Muon energy;Energy [GeV];Counts",
        bins,
        &bin_edges[0]
    );
    h_muon_energy->SetDirectory(nullptr);

#pragma omp parallel
    {
        TFile f(filename, "READ");
        TTree* t = (TTree*)f.Get("events");

        std::vector<float>* px = nullptr;
        std::vector<float>* py = nullptr;
        std::vector<float>* pz = nullptr;
        std::vector<short>* particle_type = nullptr;

        t->SetBranchAddress("px", &px);
        t->SetBranchAddress("py", &py);
        t->SetBranchAddress("pz", &pz);
        t->SetBranchAddress("particle_type", &particle_type);

        auto h_private = std::make_unique<TH1F>(
            Form("h_private_%d", omp_get_thread_num()),
            "",
            bins,
            &bin_edges[0]
        );
        h_private->SetDirectory(nullptr);

#pragma omp for
        for (int i = 0; i < nEvents; i++)
        {
            t->GetEntry(i);

            int n = px->size();

            for (int j = 0; j < n; j++)
            {
                short pid = particle_type->at(j);

                if (pid == 5 || pid == 6)
                {
                    float e = std::sqrt(
                        px->at(j) * px->at(j) +
                        py->at(j) * py->at(j) +
                        pz->at(j) * pz->at(j)
                    );

                    h_private->Fill(e);
                }
            }
        }

#pragma omp critical
        {
            h_muon_energy->Add(h_private.get());
        }

        f.Close();
    }

    // Save histogram
    std::string input_name = argv[1];
    std::string base = input_name.substr(0, input_name.find(".root"));
    
    std::string root_out = "./output/muon_histogram_" + base + ".root";
    std::string png_out  = "./output/muon_histogram_" + base + ".png";

    TFile out(root_out.c_str(), "RECREATE");
    h_muon_energy->Write();
    out.Close();

    TCanvas c("c", "Muon Energy Histogram", 900, 700);
    c.SetLogx();
    c.SetLogy();

    h_muon_energy->SetLineWidth(2);
    h_muon_energy->Draw("HIST");

    c.SaveAs(png_out.c_str());
    std::cout << "Histogram saved\n";

    double t1 = omp_get_wtime();
    std::cout << "Total runtime: " << t1 - t0 << " s\n";

    return 0;
}
