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
#include <TLine.h>
#include "TMatrixD.h"

void Sort_DSSD_PAD(TChain* chain, TChain* chainBKG, TFile* outputFile, 
    const std::string filename, const std::string isotope, const std::string Calib_filename_VMMR8, const std::string Calib_filename_MDPP32, const std::string Mapping_filename_DSSD ,const bool Substract_BKG, const std::string BKG_Calib_filename
    ){

    //Numero de modulos MDPP16, MDPP32 Y VMMR8 en nuestro experimento, en nuestro caso tenemos 1 MDPP16 con los GAGGs, dos MDPP32 donde se encuentran los PADS y el beam dump, y un VMMR8 con los DSSD. Seguramente hay más módulos y cosas conectadas a cada módulo, pero de momento no afectan a mi trabajo.
    
    //const int mdpp16_counter = 1;
    
    const int mdpp32_counter = 2;
    
    const int vmmr8_counter = 1;
    
    //Lo siguiente es incluir todos los buss que vamos a utilizar en el módulo VMMR8, que son 6
    
    const Int_t NumberOfBuss_vmmr8 = 6;

    //Numero de canales de cada módulo. No todos los canales contienen información
    
    //const Int_t NumberOfChan_mdpp16 = 16;
    
    const Int_t NumberOfChan_vmmr8 = 64;
    
    const Int_t NumberOfChan_mdpp32 = 32; 

    //Maxima multiplicidad de los eventos (cuantos canales golpean)
    
    //const Int_t maxMul_mdpp16 = 16;

    const Int_t maxMul_vmmr8 = 384; //(64*6, por si acaso tenemos un evento que da a todo lo conectado al vmmr8, que son 6 módulos,tendríamos una multiplicidad de 384, hay que estar preparado por si acaso.)
    
    const Int_t maxMul_mdpp32 = 32;
    

    TRandom3 RandomR3;

    //Definimos variables para el VMMR8
    
    Int_t mul_VMMR8[vmmr8_counter]; 
    // La multiplicidad es la multiplicidad del número de eventos totales en todo el VMMR8, por tanto depende de vmmr8_counter que es el número de VMMR8 totales que va a ser siempre 1, por lo tanto será un escalar.
    
    Int_t buss_VMMR8[vmmr8_counter][maxMul_vmmr8]; 
    int detCh_VMMR8[vmmr8_counter][maxMul_vmmr8];
    // Bus y Canal lo agrupamos ahora dependiendo de la multiplicidad. Osea, si la multiplicidad es 3, yo voy a tener 3 sucesos,
    // pues asocio a cada uno de esos sucesos su bus y su canal.

    Long64_t EventTimeStamp_VMMR8[vmmr8_counter], ExtTimeStamp_VMMR8[vmmr8_counter];
    // En un mismo entry, todos los sucesos van a tener el mismo timestamp, están sincronizados todos, por eso no depende de la multiplicidad.
    // Con el Timestamp lo que hacemos es marcar el tiempo en el que las señales de todos los módulos se leen y estas se leen a la vez.

    Long64_t values_VMMR8[vmmr8_counter][maxMul_vmmr8];
    // Lo mismo que bus y canal
    
    std::cout << "-----------------------------"  << std::endl;
    std::cout << "Reading TTree for VMMR8" << std::endl;

    for (int i = 0; i < vmmr8_counter; i++) {
        // Crea un nombre base para todas las branches
        TString baseName = TString::Format("VMMR8_%d_", i);
        chain->SetBranchAddress(baseName + "Multiplicity", &mul_VMMR8[i]);
        chain->SetBranchAddress(baseName + "Channels", detCh_VMMR8[i]);
        chain->SetBranchAddress(baseName + "Values", values_VMMR8[i]);
        chain->SetBranchAddress(baseName + "time_diff_buss", buss_VMMR8[i]);
        chain->SetBranchAddress(baseName + "Event_timestamps", &EventTimeStamp_VMMR8[i]);
        chain->SetBranchAddress(baseName + "Extended_timestamps", &ExtTimeStamp_VMMR8[i]);
    }

    // Definimos variables para el MDPP32
    
    Int_t mul_MDPP32[mdpp32_counter];
    int detCh_MDPP32[mdpp32_counter][maxMul_mdpp32];
    Long64_t values_MDPP32[mdpp32_counter][maxMul_mdpp32], EventTimeStamp_MDPP32[mdpp32_counter],ExtTimeStamp_MDPP32[mdpp32_counter], time_diff_values_MDPP32[mdpp32_counter][maxMul_mdpp32];

    // La misma logica de antes, pero ahora hay 2 mdpp32_scp, (como 2 buses), aunque tu solo vas a trabajar con 1 (de momento),
    // es mejor que tengas el código preparado para el futuro para tratar los dos.

    std::cout << "Reading TTree for MDPP32" << std::endl;
    for (int i = 0; i < mdpp32_counter; ++i) {
        TString baseName   = TString::Format("MDPP32_SCP_%d_", i);
        chain->SetBranchAddress(baseName + "Multiplicity", &mul_MDPP32[i]);
        chain->SetBranchAddress(baseName + "Channels", detCh_MDPP32[i]);
        chain->SetBranchAddress(baseName + "Values", values_MDPP32[i]);
        chain->SetBranchAddress(baseName + "Event_timestamps", &EventTimeStamp_MDPP32[i]);
        chain->SetBranchAddress(baseName + "Extended_timestamps", &ExtTimeStamp_MDPP32[i]);
    }

    //Cojo los valores de Energía
    std::vector<double> Calibration_Values = Select_Calibration_Values(isotope);
    std::vector<double> Picos_Reales = {3182.69,5156.59,5144.3,5105.5,5485.6,5442.9,5804.82,5762.7};
    std::vector<double> Ranges_values = Select_Fit_Ranges(isotope);
    const int nPeaks = Calibration_Values.size();

    //Definimos los histogramas de cada detector
    // Raw 
    char DSSD_Raw_name_Total[5][32][100]; 
    char DSSD_Raw_name_histo_Total[5][32][200]; 
    static TH1F* DSSD_Raw_hist_Total[5][32] = {nullptr}; 

    char PAD_Raw_name_Total[5][100]; 
    char PAD_Raw_name_histo_Total[5][200]; 
    static TH1F* PAD_Raw_hist_Total[5] = {nullptr};

    //Calib
    char DSSD_Calib_name_Total[5][32][100]; 
    char DSSD_Calib_name_histo_Total[5][32][200]; 
    static TH1F* DSSD_Calib_hist_Total[5][32] = {nullptr}; 

    char PAD_Calib_name_Total[5][100]; 
    char PAD_Calib_name_histo_Total[5][200]; 
    static TH1F* PAD_Calib_hist_Total[5] = {nullptr};
    
    //Check Calib
    char DSSD_Check_Calib_name[100];
    char DSSD_Check_Calib_name_histo[100];
    static TH2F* DSSD_Check_Calib_hist = {nullptr};
    
    char PAD_Check_Calib_name[100];
    char PAD_Check_Calib_name_histo[100];
    static TH2F* PAD_Check_Calib_hist = {nullptr};
    
    //Coincidencias
    char DSSD_Coinci_Two_name[5][32][100];
    char DSSD_Coinci_Two_name_histo[5][32][200];
    static TH1F* DSSD_Coinci_Two_hist[5][32] = {nullptr};
    
    //Check Calib Coincidencias
    char DSSD_Check_Calib_Coinc_name[100];
    char DSSD_Check_Calib_Coinc_name_histo[100];
    static TH2F* DSSD_Check_Calib_Coinc_hist = {nullptr};
    
    //Matriz de píxeles
    char DSSD_Coinci_Two_Pixels_name[5][100];
    char DSSD_Coinci_Two_Pixels_name_histo[5][100];
    static TH2F* DSSD_Coinci_Two_Pixels_hist[5] = {nullptr};
    
    //Diferencias entre la energía en el lado p y n
    char DSSD_Diff_EpEn_name[5][32][100];
    char DSSD_Diff_EpEn_name_histo[5][32][100];
    static TH1F* DSSD_Diff_EpEn_hist[5][32]= {nullptr};
    
    //Energía detectada por 2 strips en coincidencias
    char DSSD_Ep_En_name[5][32][100];
    char DSSD_Ep_En_name_histo[5][32][100];
    static TH2F* DSSD_Ep_En_hist[5][32] = {nullptr};
    
    char DSSD_Ep_En_name_Total[5][100];
    char DSSD_Ep_En_name_histo_Total[5][100];
    static TH2F* DSSD_Ep_En_hist_Total[5] = {nullptr};
    
    //Número de entradas en cada canal conectado del DSSD
    char DSSD_detCh_name[5][100];
    char DSSD_detCh_name_histo[5][100];
    static TH1F* DSSD_detCh_hist[5] = {nullptr};
    
    //Número de coincidencias entre strips
    char DSSD_Coinci_strips_p_name[5][100];
    char DSSD_Coinci_strips_p_name_histo[5][100];
    static TH2F* DSSD_Coinci_strips_p_hist[5] = {nullptr};

    char DSSD_Coinci_strips_n_name[5][100];
    char DSSD_Coinci_strips_n_name_histo[5][100];
    static TH2F* DSSD_Coinci_strips_n_hist[5] = {nullptr};
    
    //Histogramas después de realizar el mapping del DSSD
    char DSSD_Raw_name_Ord[5][32][100]; 
    char DSSD_Raw_name_histo_Ord[5][32][200]; 
    static TH1F* DSSD_Raw_hist_Ord[5][32] = {nullptr}; 

    //Calib
    char DSSD_Calib_name_Ord[5][32][100]; 
    char DSSD_Calib_name_histo_Ord[5][32][200]; 
    static TH1F* DSSD_Calib_hist_Ord[5][32] = {nullptr}; 
    
    //Check Calib
    char DSSD_Check_Calib_name_Ord[100];
    char DSSD_Check_Calib_name_histo_Ord[100];
    static TH2F* DSSD_Check_Calib_hist_Ord = {nullptr};
    
    //Coincidencias
    char DSSD_Coinci_Two_name_Ord[5][32][100];
    char DSSD_Coinci_Two_name_histo_Ord[5][32][200];
    static TH1F* DSSD_Coinci_Two_hist_Ord[5][32] = {nullptr};
    
    //Check Calib Coincidencias
    char DSSD_Check_Calib_Coinc_name_Ord[100];
    char DSSD_Check_Calib_Coinc_name_histo_Ord[100];
    static TH2F* DSSD_Check_Calib_Coinc_hist_Ord = {nullptr};
    
    //Matriz de píxeles
    char DSSD_Coinci_Two_Pixels_name_Ord[5][100];
    char DSSD_Coinci_Two_Pixels_name_histo_Ord[5][100];
    static TH2F* DSSD_Coinci_Two_Pixels_hist_Ord[5] = {nullptr};
    
    //Diferencias entre la energía en el lado p y n
    char DSSD_Diff_EpEn_name_Ord[5][32][100];
    char DSSD_Diff_EpEn_name_histo_Ord[5][32][100];
    static TH1F* DSSD_Diff_EpEn_hist_Ord[5][32]= {nullptr};
    
    //Energía detectada por 2 strips en coincidencias
    char DSSD_Ep_En_name_Ord[5][32][100];
    char DSSD_Ep_En_name_histo_Ord[5][32][100];
    static TH2F* DSSD_Ep_En_hist_Ord[5][32] = {nullptr};
    
    char DSSD_Ep_En_name_Total_Ord[5][100];
    char DSSD_Ep_En_name_histo_Total_Ord[5][100];
    static TH2F* DSSD_Ep_En_hist_Total_Ord[5] = {nullptr};
    
    //Número de entradas en cada canal conectado del DSSD
    char DSSD_detCh_name_Ord[5][100];
    char DSSD_detCh_name_histo_Ord[5][100];
    static TH1F* DSSD_detCh_hist_Ord[5] = {nullptr};
    
    //Número de coincidencias entre strips
    char DSSD_Coinci_strips_p_name_Ord[5][100];
    char DSSD_Coinci_strips_p_name_histo_Ord[5][100];
    static TH2F* DSSD_Coinci_strips_p_hist_Ord[5] = {nullptr};

    char DSSD_Coinci_strips_n_name_Ord[5][100];
    char DSSD_Coinci_strips_n_name_histo_Ord[5][100];
    static TH2F* DSSD_Coinci_strips_n_hist_Ord[5] = {nullptr};
    
    //A continuación creamos los histogramas calibrados para cada uno de los detectores 
    for (int i = 0; i < 5 ; i++){
        for (int j =0 ; j < 32 ; j++){
            
            sprintf(DSSD_Raw_name_Total[i][j],"DSSD_%i Strip_%i Raw", i+1,j+1); //Nombre del histograma en el archivo %i significa que ahi metes un entero
            sprintf(DSSD_Raw_name_histo_Total[i][j],"DSSD_%i Strip_%i Raw; Channel ; Counts", i+1,j+1); //Nombre del histogram ; eje X; eje Y
            DSSD_Raw_hist_Total[i][j] = new TH1F(DSSD_Raw_name_Total[i][j],DSSD_Raw_name_histo_Total[i][j],4096,0,4096); // Le metes los nombre, número de bines, bin inicial, bin final
            
            sprintf(DSSD_Calib_name_Total[i][j],"DSSD_%i Strip_%i Calib", i+1,j+1);
            sprintf(DSSD_Calib_name_histo_Total[i][j],"DSSD_%i Strip_%i Calib; Energia(KeV) ; Counts", i+1,j+1); 
            DSSD_Calib_hist_Total[i][j] = new TH1F(DSSD_Calib_name_Total[i][j],DSSD_Calib_name_histo_Total[i][j],3750,0,7500);
            
            sprintf(DSSD_Coinci_Two_name[i][j],"DSSD_%i Strip_%i Coinci", i+1,j+1);
            sprintf(DSSD_Coinci_Two_name_histo[i][j],"DSSD_%i Strip_%i Coinci; Energia(KeV) ; Counts", i+1,j+1); 
            DSSD_Coinci_Two_hist[i][j] = new TH1F(DSSD_Coinci_Two_name[i][j],DSSD_Coinci_Two_name_histo[i][j],3750,0,7500);
            
            //Se dibujan las líneas de los picos de energía para comprobar si está bien calibrado
            for(int k = 0; k < nPeaks; k ++){
                Draw_Calib_lines(DSSD_Calib_hist_Total[i][j], Calibration_Values[k],kRed);
                Draw_Calib_lines(DSSD_Coinci_Two_hist[i][j], Calibration_Values[k],kRed);
            }
            for(int l = 0; l < Picos_Reales.size(); l ++){
                Draw_Calib_lines(DSSD_Calib_hist_Total[i][j], Picos_Reales[l], 800, kBlue);
            }
            
            //Histograma con la diferencia de energía detectada en los canales n y p en el mismo evento
            sprintf(DSSD_Diff_EpEn_name[i][j]," Diff Ep-En DSSD_%i Strip_%i",i+1,j+1);
            sprintf(DSSD_Diff_EpEn_name_histo[i][j],"Ep-En Coinci DSSD_%i Strip_%i; Energia (keV); Counts",i+1,j+1);
            DSSD_Diff_EpEn_hist[i][j] = new TH1F(DSSD_Diff_EpEn_name[i][j],DSSD_Diff_EpEn_name_histo[i][j],500,-1000,1000);
            
            Draw_Calib_lines(DSSD_Diff_EpEn_hist[i][j], 0.0,8000);
        
            //Histograma 2D para comparar la energía del lado p con la del lado n
            sprintf(DSSD_Ep_En_name[i][j],"Ep vs En DSSD_%i Strip_%i",i+1,j+1);
            sprintf(DSSD_Ep_En_name_histo[i][j],"Ep vs En DSSD_%i Strip_%i; Ep (keV); En (keV)",i+1,j+1);
            DSSD_Ep_En_hist[i][j] = new TH2F(DSSD_Ep_En_name[i][j], DSSD_Ep_En_name_histo[i][j],750,0,7500,750,0,7500);
            
            Draw_Diagonal_2D(DSSD_Ep_En_hist[i][j],0.0,7500.0);
        }
    }
    
    for (int i = 0; i<5; i++){
    
        //Matriz con el número de cuentas por píxel
        sprintf(DSSD_Coinci_Two_Pixels_name[i],"N Cuentas Pixel DSSD_%i",i+1);
        sprintf(DSSD_Coinci_Two_Pixels_name_histo[i],"N Cuentas Pixel DSSD_%i; Canal p ; Canal n",i+1);
        DSSD_Coinci_Two_Pixels_hist[i] = new TH2F(DSSD_Coinci_Two_Pixels_name[i], DSSD_Coinci_Two_Pixels_name_histo[i], 16, -0.5, 15.5, 16, 15.5, 31.5);
    
        //Histograma 2D para comparar la energía del lado p con la del lado n
        sprintf(DSSD_Ep_En_name_Total[i],"Ep vs En DSSD_%i Total",i+1);
        sprintf(DSSD_Ep_En_name_histo_Total[i],"Ep vs En DSSD_%i Total; Ep (keV); En (keV)",i+1);
        DSSD_Ep_En_hist_Total[i] = new TH2F(DSSD_Ep_En_name_Total[i], DSSD_Ep_En_name_histo_Total[i],750,0,7500,750,0,7500);
            
        Draw_Diagonal_2D(DSSD_Ep_En_hist_Total[i],0.0,7500.0);
        
        //Histograma ver el número de cuentas por cada strip
        sprintf(DSSD_detCh_name[i],"N Canal, DSSD_%i",i+1);
        sprintf(DSSD_detCh_name_histo[i],"N Canal, DSSD_%i; Canal; Counts",i+1);
        DSSD_detCh_hist[i] = new TH1F(DSSD_detCh_name[i], DSSD_detCh_name_histo[i],32,0,32);
            
        //Histograma 2D para comparar el número de coincidencias entre strips del mismo lado
        sprintf(DSSD_Coinci_strips_n_name[i],"DSSD_%i Coincidencias n",i+1);
        sprintf(DSSD_Coinci_strips_n_name_histo[i],"DSSD_%i Coincidencias n; Ch; Ch",i+1);
        DSSD_Coinci_strips_n_hist[i] = new TH2F(DSSD_Coinci_strips_n_name[i], DSSD_Coinci_strips_n_name_histo[i],16,16,32,16,16,32);
    
        sprintf(DSSD_Coinci_strips_p_name[i],"DSSD_%i Coincidencias p",i+1);
        sprintf(DSSD_Coinci_strips_p_name_histo[i],"DSSD_%i Coincidencias p; Ch; Ch",i+1);
        DSSD_Coinci_strips_p_hist[i] = new TH2F(DSSD_Coinci_strips_p_name[i], DSSD_Coinci_strips_p_name_histo[i],16,0,16,16,0,16);
    }
    
    for (int i = 0; i < 5 ; i++){
            
        sprintf(PAD_Raw_name_Total[i],"PAD_%i Raw", i+1); //Nombre del histograma en el archivo %i significa que ahi metes un entero
        sprintf(PAD_Raw_name_histo_Total[i],"PAD_%i Raw; Channel ; Counts", i+1); //Nombre del histogram ; eje X; eje Y
        PAD_Raw_hist_Total[i] = new TH1F(PAD_Raw_name_Total[i],PAD_Raw_name_histo_Total[i],65536,0,65536); // Le metes los nombre, número de bines, bin inicial, bin final 
        
        sprintf(PAD_Calib_name_Total[i],"PAD_%i Calib", i+1); //Nombre del histograma en el archivo %i significa que ahi metes un entero
        sprintf(PAD_Calib_name_histo_Total[i],"PAD_%i Calib; Energia(KeV) ; Counts", i+1); //Nombre del histogram ; eje X; eje Y
        PAD_Calib_hist_Total[i] = new TH1F(PAD_Calib_name_Total[i],PAD_Calib_name_histo_Total[i],7500,0,7500); // Le metes los nombre, número de bines, bin inicial, bin final
        
        for(int k = 0; k < nPeaks; k ++){
            Draw_Calib_lines(PAD_Calib_hist_Total[i], Calibration_Values[k]);
        }
        for(int l = 0; l < Picos_Reales.size(); l ++){
            Draw_Calib_lines(PAD_Calib_hist_Total[i], Picos_Reales[l]);
        }
    
    }
    
    //Check Calib DSSD. Realizará un histograma de 2D en los cuales se verá si los picos coinciden con el valor de calibración
    sprintf(DSSD_Check_Calib_name,"Check Calib DSSD");
    sprintf(DSSD_Check_Calib_name_histo,"Check Calib DSSD; Strip ; Energia (keV)");
    DSSD_Check_Calib_hist = new TH2F(DSSD_Check_Calib_name, DSSD_Check_Calib_name_histo, 32, 0, 31,3750 , 0, 7500);
    
    //Check Calib PAD
    sprintf(PAD_Check_Calib_name,"Check Calib PAD");
    sprintf(PAD_Check_Calib_name_histo,"Check Calib PAD; Canal ; Energia (keV)");
    PAD_Check_Calib_hist = new TH2F(PAD_Check_Calib_name, PAD_Check_Calib_name_histo, 1, 5, 6, 2500, 0, 7500);
    
    //Check Calib Coinc
    sprintf(DSSD_Check_Calib_Coinc_name,"Check Calib Coinc DSSD");
    sprintf(DSSD_Check_Calib_Coinc_name_histo,"Check Calib Coinc DSSD; Strip ; Energia (keV)");
    DSSD_Check_Calib_Coinc_hist = new TH2F(DSSD_Check_Calib_Coinc_name, DSSD_Check_Calib_Coinc_name_histo, 32, 0, 31,3750, 0, 7500);
    
    
    
    for(int j = 0; j < nPeaks; j ++){
        Draw_Check_Calib_lines(DSSD_Check_Calib_hist, Calibration_Values[j]);
        Draw_Check_Calib_lines(PAD_Check_Calib_hist, Calibration_Values[j]);
        Draw_Check_Calib_lines(DSSD_Check_Calib_Coinc_hist, Calibration_Values[j]);
    }
    
        //A continuación creamos los histogramas calibrados y ordenados para cada uno de los detectores 
    for (int i = 0; i < 5 ; i++){
        for (int j =0 ; j < 32 ; j++){
            
            sprintf(DSSD_Raw_name_Ord[i][j],"DSSD_%i Strip_%i Raw Ord", i+1,j+1); //Nombre del histograma en el archivo %i significa que ahi metes un entero
            sprintf(DSSD_Raw_name_histo_Ord[i][j],"DSSD_%i Strip_%i Raw Ord; Channel ; Counts", i+1,j+1); //Nombre del histogram ; eje X; eje Y
            DSSD_Raw_hist_Ord[i][j] = new TH1F(DSSD_Raw_name_Ord[i][j],DSSD_Raw_name_histo_Ord[i][j],4096,0,4096); // Le metes los nombre, número de bines, bin inicial, bin final
            
            sprintf(DSSD_Calib_name_Ord[i][j],"DSSD_%i Strip_%i Calib Ord", i+1,j+1);
            sprintf(DSSD_Calib_name_histo_Ord[i][j],"DSSD_%i Strip_%i Calib Ord; Energia(KeV) ; Counts", i+1,j+1); 
            DSSD_Calib_hist_Ord[i][j] = new TH1F(DSSD_Calib_name_Ord[i][j],DSSD_Calib_name_histo_Ord[i][j],3750,0,7500);
            
            sprintf(DSSD_Coinci_Two_name_Ord[i][j],"DSSD_%i Strip_%i Coinci Ord", i+1,j+1);
            sprintf(DSSD_Coinci_Two_name_histo_Ord[i][j],"DSSD_%i Strip_%i Coinci Ord; Energia(KeV) ; Counts", i+1,j+1); 
            DSSD_Coinci_Two_hist_Ord[i][j] = new TH1F(DSSD_Coinci_Two_name_Ord[i][j],DSSD_Coinci_Two_name_histo_Ord[i][j],3750,0,7500);
            
            //Se dibujan las líneas de los picos de energía para comprobar si está bien calibrado
            for(int k = 0; k < nPeaks; k ++){
                Draw_Calib_lines(DSSD_Calib_hist_Ord[i][j], Calibration_Values[k],kRed);
                Draw_Calib_lines(DSSD_Coinci_Two_hist_Ord[i][j], Calibration_Values[k],kRed);
            }
            for(int l = 0; l < Picos_Reales.size(); l ++){
                Draw_Calib_lines(DSSD_Calib_hist_Ord[i][j], Picos_Reales[l], 800, kBlue);
            }
            
            //Histograma con la diferencia de energía detectada en los canales n y p en el mismo evento
            sprintf(DSSD_Diff_EpEn_name_Ord[i][j]," Diff Ep-En DSSD_%i Strip_%i Ord",i+1,j+1);
            sprintf(DSSD_Diff_EpEn_name_histo_Ord[i][j],"Ep-En Coinci DSSD_%i Strip_%i Ord; Energia (keV); Counts",i+1,j+1);
            DSSD_Diff_EpEn_hist_Ord[i][j] = new TH1F(DSSD_Diff_EpEn_name_Ord[i][j],DSSD_Diff_EpEn_name_histo_Ord[i][j],500,-1000,1000);
            
            Draw_Calib_lines(DSSD_Diff_EpEn_hist_Ord[i][j], 0.0,8000);
        
            //Histograma 2D para comparar la energía del lado p con la del lado n
            sprintf(DSSD_Ep_En_name_Ord[i][j],"Ep vs En DSSD_%i Strip_%i Ord",i+1,j+1);
            sprintf(DSSD_Ep_En_name_histo_Ord[i][j],"Ep vs En DSSD_%i Strip_%i Ord; Ep (keV); En (keV)",i+1,j+1);
            DSSD_Ep_En_hist_Ord[i][j] = new TH2F(DSSD_Ep_En_name_Ord[i][j], DSSD_Ep_En_name_histo_Ord[i][j],750,0,7500,750,0,7500);
            
            Draw_Diagonal_2D(DSSD_Ep_En_hist_Ord[i][j],0.0,7500.0);
        }
    }
    
    for (int i = 0; i<5; i++){
    
        //Matriz con el número de cuentas por píxel
        sprintf(DSSD_Coinci_Two_Pixels_name_Ord[i],"N Cuentas Pixel DSSD_%i Ord",i+1);
        sprintf(DSSD_Coinci_Two_Pixels_name_histo_Ord[i],"N Cuentas Pixel DSSD_%i Ord; Canal p ; Canal n",i+1);
        DSSD_Coinci_Two_Pixels_hist_Ord[i] = new TH2F(DSSD_Coinci_Two_Pixels_name_Ord[i], DSSD_Coinci_Two_Pixels_name_histo_Ord[i], 16, -0.5, 15.5, 16, 15.5, 31.5);
    
        //Histograma 2D para comparar la energía del lado p con la del lado n
        sprintf(DSSD_Ep_En_name_Total_Ord[i],"Ep vs En DSSD_%i Total Ord",i+1);
        sprintf(DSSD_Ep_En_name_histo_Total_Ord[i],"Ep vs En DSSD_%i Total Ord; Ep (keV); En (keV)",i+1);
        DSSD_Ep_En_hist_Total_Ord[i] = new TH2F(DSSD_Ep_En_name_Total_Ord[i], DSSD_Ep_En_name_histo_Total_Ord[i],750,0,7500,750,0,7500);
            
        Draw_Diagonal_2D(DSSD_Ep_En_hist_Total_Ord[i],0.0,7500.0);
        
        //Histograma ver el número de cuentas por cada strip
        sprintf(DSSD_detCh_name_Ord[i],"N Canal, DSSD_%i Ord",i+1);
        sprintf(DSSD_detCh_name_histo_Ord[i],"N Canal, DSSD_%i Ord; Canal; Counts",i+1);
        DSSD_detCh_hist_Ord[i] = new TH1F(DSSD_detCh_name_Ord[i], DSSD_detCh_name_histo_Ord[i],32,0,32);
            
        //Histograma 2D para comparar el número de coincidencias entre strips del mismo lado
        sprintf(DSSD_Coinci_strips_n_name_Ord[i],"DSSD_%i Coincidencias n Ord",i+1);
        sprintf(DSSD_Coinci_strips_n_name_histo_Ord[i],"DSSD_%i Coincidencias n Ord; Ch; Ch",i+1);
        DSSD_Coinci_strips_n_hist_Ord[i] = new TH2F(DSSD_Coinci_strips_n_name_Ord[i], DSSD_Coinci_strips_n_name_histo_Ord[i],16,16,32,16,16,32);
    
        sprintf(DSSD_Coinci_strips_p_name_Ord[i],"DSSD_%i Coincidencias p Ord",i+1);
        sprintf(DSSD_Coinci_strips_p_name_histo_Ord[i],"DSSD_%i Coincidencias p Ord; Ch; Ch",i+1);
        DSSD_Coinci_strips_p_hist_Ord[i] = new TH2F(DSSD_Coinci_strips_p_name_Ord[i], DSSD_Coinci_strips_p_name_histo_Ord[i],16,0,16,16,0,16);
    }
    
    //Check Calib DSSD. Realizará un histograma de 2D en los cuales se verá si los picos coinciden con el valor de calibración
    sprintf(DSSD_Check_Calib_name_Ord,"Check Calib DSSD Ord");
    sprintf(DSSD_Check_Calib_name_histo_Ord,"Check Calib DSSD Ord; Strip ; Energia (keV)");
    DSSD_Check_Calib_hist_Ord = new TH2F(DSSD_Check_Calib_name_Ord, DSSD_Check_Calib_name_histo_Ord, 32, 0, 31,3750 , 0, 7500);
    
    //Check Calib Coinc
    sprintf(DSSD_Check_Calib_Coinc_name_Ord,"Check Calib Coinc DSSD Ord");
    sprintf(DSSD_Check_Calib_Coinc_name_histo_Ord,"Check Calib Coinc DSSD Ord; Strip ; Energia (keV)");
    DSSD_Check_Calib_Coinc_hist_Ord = new TH2F(DSSD_Check_Calib_Coinc_name_Ord, DSSD_Check_Calib_Coinc_name_histo_Ord, 32, 0, 31, 3750 , 0, 7500);
    
    for(int j = 0; j < nPeaks; j ++){
        Draw_Check_Calib_lines(DSSD_Check_Calib_hist_Ord, Calibration_Values[j]);
        Draw_Check_Calib_lines(DSSD_Check_Calib_Coinc_hist_Ord, Calibration_Values[j]);
    }
    
    
    //Creamos las variables en las que se van a guardar los parámetros de calibración y guardamos estos parámetros

    double a_VMMR8[5][32], b_VMMR8[5][32];
    double a_MDPP32[5], b_MDPP32[5];
    int Ch_fis[5][32];
    
    LoadVMMR8CalibParameters(Calib_filename_VMMR8, a_VMMR8, b_VMMR8);
    
    LoadMDPP32CalibParameters(Calib_filename_MDPP32,a_MDPP32, b_MDPP32);
    
    LoadDSSDMappingParameters(Mapping_filename_DSSD,Ch_fis);

    //Cogemos el número de entradas para hacernos una idea
    Long64_t nEntries = chain->GetEntries();
    cout << "Number of Entries: " << nEntries << endl;
    int printFrequency = (nEntries > 100) ? (nEntries / 100) : 1;   // en vez de nEntries / 100 a secas

    int Loops_Time_Stamp = 0; //Por si no tengo Extended Timestamp. El TimeStamp vale para calcular el tiempo del run
    int Last_Timestamp = 0;

    //Bucle sobre las entradas
    
    int mul_DSSD5 = 0; //De momento sólo queremos la multiplicidad de los eventos que han dado en el DSSD 5
    int n_DSSD5 = 0; //Necesitamos un n y un p por separado, por lo que creamos estas variables para comprobar
    int p_DSSD5 = 0;
    int k_p = 0; //Definimos estos enteros en los que se va a guardar el número de multiplicidad que corresponde a cada coincidencia
    int k_n = 0;
    double E_p;
    double E_n;
    std::vector <int> idx_p;
    std::vector <int> idx_n;
    
    for (int i = 0; i < int(nEntries); i++){

        double Energy_DSSD[5][32] = {0}; // En teoría los detectores solo pueden coger una señal por evento, por eso es un escalar
        double Energy_PAD[5] = {0};
        
        //Variables para el histograma de coincidencias
        mul_DSSD5 = 0; //Inicializamos la multiplicidad en el detector 5
        n_DSSD5 = 0;
        p_DSSD5 = 0;
        idx_p.clear();
        idx_n.clear();
        
        chain -> GetEvent(i);
        
        if(Last_Timestamp > EventTimeStamp_MDPP32[0]){ // El TimeStamp debería ser el mismo para todos los módulos, así que con coger uno vale
            Loops_Time_Stamp++;
        }
        Last_Timestamp = EventTimeStamp_MDPP32[0];

        for (int j = 0; j < mdpp32_counter; j++){ //Bucle de los MDPP32

            for (Int_t k = 0; k < mul_MDPP32[j]; k++){ // mul_MDPP32[j] nos da el número de eventos registrados en el módulo para esta entrada
                int PAD_Number = detCh_MDPP32[j][k]; //Se obtiene el canal físico del evento (el detector que lo ha medido)
                if (PAD_Number == 5){ //De momento sólo voy a analizar el PAD 5
                    PAD_Raw_hist_Total[PAD_Number-1]->Fill(values_MDPP32[j][k]); //Se llena el histograma Raw

                    double CalibratedEnergy_MDPP32 = (a_MDPP32[PAD_Number-1]*(values_MDPP32[j][k]+RandomR3.Uniform()-0.5)+b_MDPP32[PAD_Number-1]); //Se aplica el ajuste lineal de la calibración

                    if(CalibratedEnergy_MDPP32 > 100.0){ //Ponemos un Noise Cut 
                        
                        //Llenamos el histograma calibrado con el parámetro obtenido de la calibración.
                        Energy_PAD[PAD_Number-1] = CalibratedEnergy_MDPP32;
                        PAD_Calib_hist_Total[PAD_Number-1]->Fill(CalibratedEnergy_MDPP32);
                        
                        PAD_Check_Calib_hist->Fill((double)PAD_Number, CalibratedEnergy_MDPP32);

                    }   
                }
            }
        }
        
        for (int j = 0; j < vmmr8_counter; j++){// Bucle por si hubiera más VMMR8
            
            for (Int_t k = 0; k < mul_VMMR8[j]; k++){
                    
                int BUSS_Number = buss_VMMR8[j][k];
                if (BUSS_Number == 5){ //Si no se encuentra en el buss al que está conectado el DSSD que nos interesa pasa al siguiente.
                    
                    int idx= detCh_VMMR8[j][k];
                    if (idx < 0 || idx > 47) continue;
                    if (idx > 15 && idx < 32) continue; //Nos quedamos sólo con los canales con strips de los DSSD conectadas
                    
                    if (idx > 15){ //Evitamos irnos fuera del rango del array (max 31)
                        idx=idx-16; 
                    }
                    DSSD_Raw_hist_Total[4][idx]->Fill(values_VMMR8[j][k]);//Values son los canales electrónicos, detCh los físicos (detectores)                                     
                    DSSD_Raw_hist_Ord[4][Ch_fis[4][idx]]->Fill(values_VMMR8[j][k]);//Values son los canales electrónicos, detCh los físicos (detectores)

                    
                    double CalibratedEnergy_VMMR8 = (a_VMMR8[4][idx]*(values_VMMR8[j][k]+RandomR3.Uniform()-0.5)+b_VMMR8[4][idx]); //Se aplica el ajuste lineal de la calibración
                    
                    if(CalibratedEnergy_VMMR8 > 350.0){ //Ponemos un Noise Cut 
                                                
                       if (buss_VMMR8[j][k] == 5 && detCh_VMMR8[j][k] >= 0 && detCh_VMMR8[j][k] <= 15){//Coindiciones para que se detecte en una strip p
                    
                            mul_DSSD5++;//Se actualizan los parámetros necesarios
                            p_DSSD5++;
                            E_p = CalibratedEnergy_VMMR8;
                            idx_p.push_back(idx);
                        }
                       
                        else if (buss_VMMR8[j][k] == 5 && detCh_VMMR8[j][k] >= 32 && detCh_VMMR8[j][k] <= 47){ //Coindiciones para que se detecte en una strip n
                           
                            mul_DSSD5++;//Si forma parte de los canales que pertenecen al detector se aumenta su multiplicidad
                            n_DSSD5++;
                            E_n=CalibratedEnergy_VMMR8;
                            idx_n.push_back(idx);
                       
                        }
                        //Llenamos el histograma calibrado con el parámetro obtenido de la calibración.
                        Energy_DSSD[4][idx] = CalibratedEnergy_VMMR8;
                        DSSD_Calib_hist_Total[4][idx]->Fill(CalibratedEnergy_VMMR8);
                        DSSD_Calib_hist_Ord[4][Ch_fis[4][idx]]->Fill(CalibratedEnergy_VMMR8);

                        
                        DSSD_Check_Calib_hist->Fill((double)idx, CalibratedEnergy_VMMR8);
                        DSSD_Check_Calib_hist_Ord->Fill((double)Ch_fis[4][idx], CalibratedEnergy_VMMR8);
                    } 
                }
            }
        }
        
        if ( mul_DSSD5 == 2 && p_DSSD5*n_DSSD5 == 1 ){//Se comprueba que hay multiplicidad 2 y que hay una strip n y una p
            
            DSSD_Diff_EpEn_hist[4][idx_p[0]]->Fill(E_p-E_n);
            DSSD_Diff_EpEn_hist[4][idx_n[0]]->Fill(E_p-E_n);
            DSSD_Ep_En_hist[4][idx_p[0]]->Fill(E_p,E_n);//Llenamos el histograma 2D para comparar la Ep y la En
            DSSD_Ep_En_hist[4][idx_n[0]]->Fill(E_p,E_n);//Llenamos el histograma 2D para comparar la Ep y la En            
            DSSD_Ep_En_hist_Total[4]->Fill(E_p,E_n);
                        
            DSSD_Diff_EpEn_hist_Ord[4][Ch_fis[4][idx_p[0]]]->Fill(E_p-E_n);
            DSSD_Diff_EpEn_hist_Ord[4][Ch_fis[4][idx_n[0]]]->Fill(E_p-E_n);
            DSSD_Ep_En_hist_Ord[4][Ch_fis[4][idx_p[0]]]->Fill(E_p,E_n);//Llenamos el histograma 2D para comparar la Ep y la En
            DSSD_Ep_En_hist_Ord[4][Ch_fis[4][idx_n[0]]]->Fill(E_p,E_n);//Llenamos el histograma 2D para comparar la Ep y la En            
            DSSD_Ep_En_hist_Total_Ord[4]->Fill(E_p,E_n);
            
            if ( std::abs(E_p-E_n) <= 200 ){ //Tienen que dar la misma energía los canales p y n (2sigmas)
                
                DSSD_Coinci_Two_hist[4][idx_p[0]]->Fill(E_p);//Llenar los histogramas p y n
                DSSD_Coinci_Two_hist[4][idx_n[0]]->Fill(E_n);
                
                DSSD_Coinci_Two_Pixels_hist[4]->Fill((double)idx_p[0],(double)idx_n[0]);//Llenamos el histograma 2D con el número de cuentas por píxel
                DSSD_Check_Calib_Coinc_hist->Fill((double)idx_p[0],E_p);
                DSSD_Check_Calib_Coinc_hist->Fill((double)idx_n[0],E_n);
                DSSD_detCh_hist[4]->Fill((double)idx_n[0]);
                DSSD_detCh_hist[4]->Fill((double)idx_p[0]);//Llenamos el histograma con el número de cuentas por cada strip

                DSSD_Coinci_Two_hist_Ord[4][Ch_fis[4][idx_p[0]]]->Fill(E_p);//Llenar los histogramas p y n
                DSSD_Coinci_Two_hist_Ord[4][Ch_fis[4][idx_n[0]]]->Fill(E_n);
                
                DSSD_Coinci_Two_Pixels_hist_Ord[4]->Fill(Ch_fis[4][idx_p[0]],Ch_fis[4][idx_n[0]]);//Llenamos el histograma 2D con el número de cuentas por píxel
                DSSD_Check_Calib_Coinc_hist_Ord->Fill(Ch_fis[4][idx_p[0]],E_p);
                DSSD_Check_Calib_Coinc_hist_Ord->Fill(Ch_fis[4][idx_n[0]],E_n);
                DSSD_detCh_hist_Ord[4]->Fill(Ch_fis[4][idx_n[0]]);
                DSSD_detCh_hist_Ord[4]->Fill(Ch_fis[4][idx_p[0]]);//Llenamos el histograma con el número de cuentas por cada strip
            }
        }
        if (mul_DSSD5 >= 2 && p_DSSD5 == 2){//Vamos a rellenar el histograma con los eventos que han dado entre dos strips p
            DSSD_Coinci_strips_p_hist[4]->Fill((double)idx_p[0],(double)idx_p[1]);
            DSSD_Coinci_strips_p_hist[4]->Fill((double)idx_p[1],(double)idx_p[0]);
 
            DSSD_Coinci_strips_p_hist_Ord[4]->Fill(Ch_fis[4][idx_p[0]],Ch_fis[4][idx_p[1]]);
            DSSD_Coinci_strips_p_hist_Ord[4]->Fill(Ch_fis[4][idx_p[1]],Ch_fis[4][idx_p[0]]);
        }

        if (mul_DSSD5 >= 2 && n_DSSD5 == 2){//Vamos a rellenar el histograma con los eventos que han dado entre dos strips n
            DSSD_Coinci_strips_n_hist[4]->Fill((double)idx_n[0],(double)idx_n[1]);
            DSSD_Coinci_strips_n_hist[4]->Fill((double)idx_n[1],(double)idx_n[0]);
            
            DSSD_Coinci_strips_n_hist_Ord[4]->Fill(Ch_fis[4][idx_n[0]],Ch_fis[4][idx_n[1]]);
            DSSD_Coinci_strips_n_hist_Ord[4]->Fill(Ch_fis[4][idx_n[1]],Ch_fis[4][idx_n[0]]);
        }
    
        if (i % printFrequency == 0 || i == nEntries - 1) {
        double progress = (double(i) / nEntries) * 100.0;
        std::cout << " Creando histos Raw, Calib y Coinci: " << progress << "% \r";
        std::cout.flush();
        }
    }
    /*Lo utilicé para estudiar por cuál de las strips p comienza
    //Vamos a ordenar las strips p utilizando el número de cuentas
    std::vector <int> Cuentas_p;
    std::vector <int> Cuentas_p_ord;
    std::vector <int> orden;
    Cuentas_p.clear();
    Cuentas_p_ord.clear();
    orden.clear();
    for (int i = 0; i < 16; i++){//Se guarda en el vector el número de cuentas en cada strip
        Cuentas_p.push_back(DSSD_detCh_hist[4]->GetBinContent(i+1));//Se 
        std::cout<<Cuentas_p[i]<<std::endl;
    }        
    Cuentas_p_ord=Cuentas_p;
    std::sort(Cuentas_p_ord.begin(),Cuentas_p_ord.end());//Ordenamos de menor a mayor el número de cuentas
    std::reverse(Cuentas_p_ord.begin(),Cuentas_p_ord.end());//Damos la vuelta al vector
    std::cout<<""<<std::endl;
    for (int i = 0; i < 16; i++){
         std::cout<<Cuentas_p_ord[i]<<std::endl;
    }
    std::cout<<""<<std::endl;
    for (int i = 0; i < 16; i++){//Creamos un vector en el que se guarde la posición de las strips
        for (int j = 0; j < 16; j++){
            if (Cuentas_p_ord[i] == Cuentas_p[j]){
                orden.push_back(j+1);
                std::cout<<orden[i]<<std::endl;
            }
        }
    }
    */

    chain->GetEvent(nEntries);//Vamos a hacer que aparezca por pantalla el tiempo que ha estado midiendo.
    std::cout << "Time Stamp: " << EventTimeStamp_MDPP32[0] << "\n";
    std::cout << "Extended Time Stamp: " << ExtTimeStamp_MDPP32[0] << "\n";
    std::cout << "Loops: " << Loops_Time_Stamp << "\n";
    double Time;
    Time = (EventTimeStamp_MDPP32[0] + Loops_Time_Stamp*pow(2,30))/16.0; // En us
    double Time_hours = Time*pow(10,-6)/3600.0;
    std::cout<<"Ha medido "<< Time_hours << " horas"<<"\n";
    
    
    // Vamos a guardar los histogramas obtenidos
    gROOT->SetBatch(kTRUE);

    TDirectory* DSSDsDir = outputFile->mkdir("DSSDs");
    Save_DSSD_Histos(DSSDsDir, "Raw", DSSD_Raw_hist_Total);
    Save_DSSD_Histos(DSSDsDir, "Calib", DSSD_Calib_hist_Total);

    TDirectory* PADsDir = outputFile->mkdir("PADs");
    Save_PAD_Histos(PADsDir, "Raw", PAD_Raw_hist_Total);
    Save_PAD_Histos(PADsDir, "Calib", PAD_Calib_hist_Total);
    
    TDirectory* SortsDir = outputFile->mkdir("Sort");
    Save_Histo(DSSD_Check_Calib_hist, SortsDir);
    Save_Histo(PAD_Check_Calib_hist, SortsDir);
    
    TDirectory* CoincisDSSDDir = outputFile->mkdir("CoincisDSSD");
    Save_DSSD_Histos(CoincisDSSDDir, "Strips", DSSD_Coinci_Two_hist);
    Save_DSSD_Histos(CoincisDSSDDir, "Diff_EpEn", DSSD_Diff_EpEn_hist);
    Save_DSSD_Histos_2D(CoincisDSSDDir, "EpvsEn", DSSD_Ep_En_hist);

    // Creamos los subdirectorios directamente desde CoincisDSSDDir
    TDirectory* PixelsDir = CoincisDSSDDir->mkdir("Pixels");
    TDirectory* CheckCalibCoincDir = CoincisDSSDDir->mkdir("CheckCalibCoinci");
    
    // Obtenemos el puntero al directorio EpvsEn creado previamente por Save_DSSD_Histos_2D
    TDirectory* EpvsEnDir = CoincisDSSDDir->GetDirectory("EpvsEn");

    if (CheckCalibCoincDir) Save_Histo(DSSD_Check_Calib_Coinc_hist, CheckCalibCoincDir);

    for (int i = 0; i < 5; i++) {
        if (PixelsDir) Save_Histo(DSSD_Coinci_Two_Pixels_hist[i], PixelsDir);
        if (EpvsEnDir) Save_Histo(DSSD_Ep_En_hist_Total[i], EpvsEnDir);
        if (SortsDir)  Save_Histo(DSSD_detCh_hist[i], SortsDir);
        if (SortsDir) Save_Histo(DSSD_Coinci_strips_n_hist[i], SortsDir);
        if (SortsDir) Save_Histo(DSSD_Coinci_strips_p_hist[i], SortsDir);
    }
    
    TDirectory* DSSDOrdsDir = outputFile->mkdir("DSSDs Ord");
    Save_DSSD_Histos(DSSDOrdsDir, "Raw Ord", DSSD_Raw_hist_Ord);
    Save_DSSD_Histos(DSSDOrdsDir, "Calib Ord", DSSD_Calib_hist_Ord);
    
    TDirectory* SortOrdsDir = outputFile->mkdir("Sort Ord");
    Save_Histo(DSSD_Check_Calib_hist_Ord, SortOrdsDir);
    
    TDirectory* CoincisDSSDOrdDir = outputFile->mkdir("CoincisDSSD Ord");
    Save_DSSD_Histos(CoincisDSSDOrdDir, "Strips Ord", DSSD_Coinci_Two_hist_Ord);
    Save_DSSD_Histos(CoincisDSSDOrdDir, "Diff_EpEn Ord", DSSD_Diff_EpEn_hist_Ord);
    Save_DSSD_Histos_2D(CoincisDSSDOrdDir, "EpvsEn Ord", DSSD_Ep_En_hist_Ord);

    // Creamos los subdirectorios directamente desde CoincisDSSDDir
    TDirectory* PixelOrdsDir = CoincisDSSDOrdDir->mkdir("Pixels Ord");
    TDirectory* CheckCalibCoincOrdDir = CoincisDSSDOrdDir->mkdir("CheckCalibCoinci Ord");
    
    // Obtenemos el puntero al directorio EpvsEn creado previamente por Save_DSSD_Histos_2D
    TDirectory* EpvsEnOrdDir = CoincisDSSDOrdDir->GetDirectory("EpvsEn Ord");

    if (CheckCalibCoincOrdDir) Save_Histo(DSSD_Check_Calib_Coinc_hist_Ord, CheckCalibCoincOrdDir);

    for (int i = 0; i < 5; i++) {
        if (PixelOrdsDir) Save_Histo(DSSD_Coinci_Two_Pixels_hist_Ord[i], PixelOrdsDir);
        if (EpvsEnOrdDir) Save_Histo(DSSD_Ep_En_hist_Total_Ord[i], EpvsEnOrdDir);
        if (SortOrdsDir)  Save_Histo(DSSD_detCh_hist_Ord[i], SortOrdsDir);
        if (SortOrdsDir) Save_Histo(DSSD_Coinci_strips_n_hist_Ord[i], SortOrdsDir);
        if (SortOrdsDir) Save_Histo(DSSD_Coinci_strips_p_hist_Ord[i], SortOrdsDir);
    }
}
