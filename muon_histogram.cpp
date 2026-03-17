#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

#include <vector>
#include <iostream>
#include <cmath>
#include <omp.h>

#include <TCanvas.h>
#include <TStyle.h>
#include <memory>
#include <TString.h>

#include <TClass.h>
#include <TParameter.h>

int main(int argc, char* argv[]) {

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <root_file>\n";
        return 1;
    }

    double t0 = omp_get_wtime();

    // these lines are necessary to prevent multithreading crush when reading root files
    std::string filepath = std::string("./data/") + argv[1];
    const char* filename = filepath.c_str();
    ROOT::EnableThreadSafety();

    gROOT->ProcessLine("#include <vector>");

    TFile warm_file(filename, "READ");
    TTree* warm_tree = (TTree*)warm_file.Get("events");

    std::vector<float>* warm_energy = nullptr;
    std::vector<short>* warm_particle_type = nullptr;

    warm_tree->SetBranchAddress("energy", &warm_energy);
    warm_tree->SetBranchAddress("particle_type", &warm_particle_type);
    warm_tree->GetEntry(0);

    warm_file.Close();

    // after initializing root to multithreading well-functioning we open our file
    TFile file(filename, "READ");
    TTree* tree = (TTree*)file.Get("events");

    int nEvents = tree->GetEntries();

    // --- Read primary info
    auto primary_energy  = (TParameter<float>*)file.Get("primary_energy");
    auto primary_zenith  = (TParameter<float>*)file.Get("primary_zenith");
    auto primary_azimuth = (TParameter<float>*)file.Get("primary_azimuth");
    auto primary_type    = (TParameter<int>*)file.Get("primary_type");

    std::cout << "=== PRIMARY INFO ===\n";
    if (primary_energy)
        std::cout << "Energy  = " << primary_energy->GetVal() << std::endl;
    if (primary_zenith)
        std::cout << "Zenith  = " << primary_zenith->GetVal() << std::endl;
    if (primary_azimuth)
        std::cout << "Azimuth = " << primary_azimuth->GetVal() << std::endl;
    if (primary_type)
        std::cout << "Type    = " << primary_type->GetVal() << std::endl;
    
    return 0;
    //find energy range (for histogram bin) and total muons
    double global_min = 1e30;
    double global_max = -1e30;
    long long total_muons = 0;

// multithreading to load energy and particle_type from all events    
#pragma omp parallel
{
    TFile f(filename);
    TTree* t = (TTree*)f.Get("events");

    std::vector<float>* energy = nullptr;
    std::vector<short>* particle_type = nullptr;

    t->SetBranchAddress("energy", &energy);
    t->SetBranchAddress("particle_type", &particle_type);

    double local_min = 1e30;
    double local_max = -1e30;
    long long local_count = 0;

// event paralellization 
#pragma omp for schedule(dynamic,1)
    for(int i = 0; i < nEvents; i++)
    {
        t->GetEntry(i);

        int n = energy->size();

        for(int j = 0; j < n; j++)
        {
            short pid = particle_type->at(j);

            if(pid == 5 || pid == 6)
            {
                float e = energy->at(j);

                if(e < local_min) local_min = e;
                if(e > local_max) local_max = e;

                local_count++;
            }
        }
    }

// update global energy min/max 
#pragma omp critical
    {
        if(local_min < global_min) global_min = local_min;
        if(local_max > global_max) global_max = local_max;
        total_muons += local_count;
    }
}

    std::cout << "Total muons = " << total_muons << std::endl;

    // --- Compute number of bins

    int bins = static_cast<int>(sqrt(total_muons));

    std::cout << "Histogram bins = " << bins << std::endl;

    // --- Logarithmic bin edges

    std::vector<double> bin_edges(bins + 1);

    double log_min = log10(global_min);
    double log_max = log10(global_max);

    for(int i = 0; i <= bins; i++)
    {
        bin_edges[i] = pow(10, log_min + (log_max - log_min) * i / bins);
    }

    // create histogram

    auto h_muon_energy = std::make_unique<TH1F>(
    "h_muon_energy",
    "Muon energy;Energy [GeV];Counts",
    bins,
    &bin_edges[0]
    );
    h_muon_energy->SetDirectory(nullptr);

   // fill histogram

#pragma omp parallel
{
    TFile f(filename);
    TTree* t = (TTree*)f.Get("events");

    std::vector<float>* energy = nullptr;
    std::vector<short>* particle_type = nullptr;

    t->SetBranchAddress("energy", &energy);
    t->SetBranchAddress("particle_type", &particle_type);

    auto h_private = std::make_unique<TH1F>(
        Form("h_private_%d", omp_get_thread_num()),
        "",
        bins,
        &bin_edges[0]
    );
    h_private->SetDirectory(nullptr);

#pragma omp for
    for(int i = 0; i < nEvents; i++)
    {
        t->GetEntry(i);

        int n = energy->size();

        for(int j = 0; j < n; j++)
        {
            short pid = particle_type->at(j);

            if(pid == 5 || pid == 6)
            {
                h_private->Fill(energy->at(j));
            }
        }
    }

#pragma omp critical
    {
        h_muon_energy->Add(h_private.get());
    }

    f.Close();
}

    // --- Save histogram
    
    std::string root_out = Form("./output/muon_histogram_%s", argv[1]);
    std::string png_out  = Form("./output/muon_histogram_%s.png", argv[1]);
    
    TFile out(root_out.c_str(),"RECREATE");
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

