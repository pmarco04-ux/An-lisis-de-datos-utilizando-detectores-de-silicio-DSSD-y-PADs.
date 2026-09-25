// Empezamos metiendo las librerias necesarias 
#include <string>
#include <utility> 

#include <fstream>
#include <sstream>

#include <string>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <filesystem> 

// ROOT libraries
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom3.h>
#include <TLatex.h>
#include <TSpectrum.h>
#include <TVirtualFitter.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TTimer.h>
#include <TROOT.h>
#include <TStopwatch.h>
#include "TMatrixD.h"

// Incluir códigos extra
#include "General_DSSD_PAD.cpp"
#include "Calibracion_General_DSSD_PAD.cpp"
#include "Sort_DSSD_PAD.cpp"

int Cal_Sort;
bool Calibrar, Sort;


std::cout<<"¿Quieres calibrar o realizar el análisis?"<<std::endl<< "1. Calibrar"<<std::endl<<std::endl<<"2. Sort"<<std::endl; //Pregunta qué tipo de detector quieres calibrar para agilizar el proceso y no tener que calibrar todo cada vez.
            
std::cin>>Cal_Sort;

if (Cal_Sort == 1){
    Calibrar = 1;
    Sort = 0;
}

else {
    Calibrar = 0;
    Sort = 1;
}

int DSSD_PAD_Analysis(){ 
    gStyle->SetOptStat(1111111); 
    
    std::vector<std::string> Runs_Analizar = {"/home/daniferu/Desktop/Master_Pablo/Root_Alphas/After/IS690B_Calb_run006_251021_195318_part001_4.root"};
    
    
    
    int Amount_of_Runs = Runs_Analizar.size();

    for(int i = 0; i < Amount_of_Runs; i++){
        TChain* eventChain = new TChain("EventTree");

        //Cojo la ruta de los archivos a analizar
        std::string direction = Runs_Analizar[i]; 
        size_t rootPos = direction.find(".root"); 
        size_t backslashPos = direction.rfind('/', rootPos);
        std::string Name = direction.substr(backslashPos + 1, rootPos - backslashPos - 1);
        eventChain->Add(direction.c_str());
        std::cout << "Procesing "<< Name << ".root" << "\n";
        std::cout << "-------------------------- "<< "\n";

        if(Calibrar){
            std::vector<std::string> Calib_Nombre ={"_Calib.root"};
            std::string Appendix_Calib = Calib_Nombre[0];
            std::string outputFileName_Calib = "./Analized/"+ Name + Appendix_Calib; 
            TFile* outputFile_Calib = new TFile(outputFileName_Calib.c_str(),"recreate");
            std::cout << "Output Name: "<< outputFileName_Calib  << "\n";
            std::cout << "-------------------------- " << "\n"; 

            TStopwatch timer;
            timer.Start();
            
            std::string isotope;
            int CalibDetector = 0;
            int numCalibDetector = 0;
    
            std::cout<<"¿Qué tipo de detector quieres calibrar? Introduce el número de la opción deseada"<<std::endl<< "0. Todos"<<std::endl<<std::endl<<"1. GAGG"<<std::endl<<"2. DSSD"<<std::endl<<"3. PAD"<<std::endl; //Pregunta qué tipo de detector quieres calibrar para agilizar el proceso y no tener que calibrar todo cada vez.
            
            std::cin>>CalibDetector;
            
            if (CalibDetector == 1){
                std::cout<<"¿Cuál de todos los GAGGs quieres calibrar? Introduce la opción deseada"<<std::endl<<"0. (Todos)" <<std::endl<<"1 "<<std::endl<<"2 "<<std::endl<<"3 "<<std::endl<<"4 "<<std::endl<<"5 "<<std::endl<<"6 "<<std::endl<<"7 "<<std::endl<<"8 ";
                std::cin>>numCalibDetector;
                
                if (numCalibDetector > 8){
                    numCalibDetector = 0;
                    std::cout<<"Valor introducido no válido, se calibrarán todos los GAGGs"<<std::endl;
                }
            }
                
            else if (CalibDetector == 2){
                std::cout<<"¿Cuál de todos los DSSD quieres calibrar? Introduce la opción deseada" <<std::endl<< "0. (Todos)"<<std::endl<<"1 "<<std::endl<<"2 "<<std::endl<<"3 "<<std::endl<<"4 "<<std::endl<<"5 ";
                std::cin>>numCalibDetector;
                
                if (numCalibDetector > 5){
                    numCalibDetector = 0;
                    std::cout<<"Valor introducido no válido, se calibrarán todos los DSSD"<<std::endl;
                }
                std::string N_run;
                std::cout<<"¿Qué archivo de calibración quieres utilizar?"<<std::endl<<"1. part001"<<std::endl<<"2. part002"<<std::endl<<"3. part003"<<std::endl<<"4. part004"<<std::endl<<"1_4. part001_4"<<std::endl;
                std::cin>>N_run;
                Runs_Analizar = {"/home/daniferu/Desktop/Master_Pablo/Root_Alphas/After/IS690B_Calb_run006_251021_195318_part00"+N_run+".root"};
                
                isotope = "Alpha DSSD";

            }
                
            else if (CalibDetector == 3){
                std::cout<<"¿Cuál de todos los PAD quieres calibrar? Introduce la opción deseada"<<std::endl<< "0. (Todos)" <<std::endl<<"1 "<<std::endl<<"2 "<<std::endl<<"3 "<<std::endl<<"4 "<<std::endl<<"5 ";
                std::cin>>numCalibDetector;
                if (numCalibDetector > 5){
                    numCalibDetector = 0;
                    std::cout<<"Valor introducido no válido, se calibrarán todos los PAD"<<std::endl;
                }
                Runs_Analizar = { 
                "/home/daniferu/Desktop/Master_Pablo/Root_Alphas/After/IS690B_Calb_run007_251022_113315_part001.root"};
                
                isotope = "Alpha PAD";

            }
            
            else {
                CalibDetector = 0;
                numCalibDetector = 0;
                std::cout<<"El valor introducido no es válido, se calibrarán todos los detectores."<<std::endl;
            }
            
            std::string NoiseCutRawFile = "Calibracion/Noise_Cut_Raw_DSSD.txt";
                            
            Calibracion_General_DSSD_PAD(eventChain, outputFile_Calib, Name, isotope, NoiseCutRawFile, CalibDetector, numCalibDetector);

            timer.Stop();
            timer.Print();
        }
    
        if(Sort){
            std::vector<std::string> Nombre ={"_Analized.root"}; 
            TChain* eventChainBKG = new TChain("EventTree");
            std::string Appendix = Nombre[0];
            std::string outputFileName_Analized = "./Analized/"+ Name + Appendix;
            TFile* outputFile_Analized = new TFile(outputFileName_Analized.c_str(),"recreate");
            std::cout << "Output Name: "<< outputFileName_Analized  << "\n";
            std::cout << "-------------------------- " << "\n";
            
            eventChainBKG->Add("ROOT/Practicas_Pablo_BKGrun001_260306_161317_part001.root");

            //bool Substract_BKG = true;
            bool Substract_BKG = false;
            
            std::string isotope = "Alpha DSSD";

            std::string Calib_filename_MDPP32 = "Calibracion/Calibration_PAD_IS690B_Calb_run007_251022_113315_part001.txt";
            std::string Calib_filename_VMMR8 = "Calibracion/Calibration_DSSD_IS690B_Calb_run006_251021_195318_part001_4.txt";
            std::string Mapping_filename_DSSD = "Calibracion/Mapping_DSSD_IS690B.txt";
            std::string BKG_Calib_filename = "Calibracion/Calibration_Practicas_Pablo_22Na_Coincidencias_SetUp1_run001_260317_112856_part001.txt";
            TStopwatch timer;
            timer.Start();

            Sort_DSSD_PAD(eventChain, eventChainBKG, outputFile_Analized, Name, isotope, Calib_filename_VMMR8, Calib_filename_MDPP32,Mapping_filename_DSSD, Substract_BKG, BKG_Calib_filename);

            timer.Stop();
            timer.Print();        
        }                                                                                                                                                                                                                                                                           
    }

    //Para poder visualizar lo reactiva porque arriba lo hemos desactivado
    gROOT->SetBatch(kFALSE);
    
    // Tiene que devolver un entero y normalmente se pone 0 porque eso significa que todo fue bien
    
    return 0;
}

