#include <string>
#include <utility>

#include <fstream>
#include <sstream>

#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <filesystem>

//Librerias de ROOT
#include <TChain.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TTree.h>
#include <TFile.h>
#include <TRandom.h>
#include <TLatex.h>
#include <TSpectrum.h>
#include <TVirtualFitter.h>
#include <TSystem.h>
#include <TLine.h>

//Veamos como leemos los datos. Esta función es un void porque no devuelve nada solo analiza/cambia los paránmetros que le metemos
void Calibracion_General_DSSD_PAD(TChain* chain,TFile* outputFile, std::string filename, std::string isotope, std::string NoiseCutRawFile, int CalibDetector, int numCalibDetector
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

    //Defino los histogramas de los DSSD, habrá que calibrar cada strip por separado, por lo que serán 32 calibraciones por cada DSSD
    
    char DSSD_Raw_name_Total[5][32][100]; //Nombre en el archivo outputFile. Son caracteres. [5] indica el número de detectores, [32] el número de strips por cada detector (5x32) histogramas y [100] la longitud de la cadena de caracteres
    char DSSD_Raw_name_histo_Total[5][32][200]; //Nombre del histograma
    static TH1F* DSSD_Raw_hist_Total[5][32]= {nullptr}; //Histogramas 1 dimensional TH1 y la F es para floats. Lo inicializo a nullptr que significa puntero nulo para evitar fallos

    //Ahora hay que darle nombre a los histogramas, hay que tener en cuenta que los números de los canales y los detectores no tienen porqué ser los mismos que los físicos. Habrá que ordenarlos en el análisis.
    for (int i = 0; i < 5 ; i++){
        for (int j =0 ; j < 32 ; j++){
            
            sprintf(DSSD_Raw_name_Total[i][j],"DSSD_%i Strip_%i Raw", i+1,j+1); //Nombre del histograma en el archivo %i significa que ahi metes un entero
            sprintf(DSSD_Raw_name_histo_Total[i][j],"DSSD_%i Strip_%i Raw; Channel ; Counts", i+1,j+1); //Nombre del histogram ; eje X; eje Y
            DSSD_Raw_hist_Total[i][j] = new TH1F(DSSD_Raw_name_Total[i][j],DSSD_Raw_name_histo_Total[i][j],4096,0,4096); // Le metes los nombre, número de bines, bin inicial, bin final   
        
        }
    }
    
    //Defino los histogramas de los PADs
    
    char PAD_Raw_name_Total[5][100]; //Nombre en el archivo outputFile. Son caracteres. [8] indica el número de histogramas y [100] la longitud de la cadena de caracteres
    char PAD_Raw_name_histo_Total[5][200]; //Nombre del histograma
    static TH1F* PAD_Raw_hist_Total[5]= {nullptr}; //Histogramas 1 dimensional TH1 y la F es para floats. Lo inicializo a nullptr que significa puntero nulo para evitar fallos

    //Ahora hay que darle nombre a los histogramas
    for (int i = 0; i < 5 ; i++){
        
        sprintf(PAD_Raw_name_Total[i],"PAD_%i Raw", i+1); //Nombre del histograma en el archivo %i significa que ahi metes un entero
        sprintf(PAD_Raw_name_histo_Total[i],"PAD_%i Raw; Channel ; Counts", i+1); //Nombre del histogram ; eje X; eje Y
        PAD_Raw_hist_Total[i] = new TH1F(PAD_Raw_name_Total[i],PAD_Raw_name_histo_Total[i],65536,0,65536); // Le metes los nombre, número de bines, bin inicial, bin final   
    
    }
    
    //Cargamos los parámetros de noise cut guardados
    int noise_cut_dig[5][32], noise_cut_fis[5][32];
    LoadDSSDNoise_CutRawParameters(NoiseCutRawFile,noise_cut_dig, noise_cut_fis);
    
    //Cogemos el número de entradas para hacernos una idea
    Long64_t nEntries = chain->GetEntries();
    std::cout << "Number of Entries: " << nEntries << std::endl;
   
    //Creamos los histogramas Raw de los DSSD y los PAD
    
    //De momento sólo me interesa el DSSD_5, por lo que me quedo sólo con el 5
    
    int nDSSD=0; //Número de histograma que corresponde al detector
    int n_buss=0; //Número de buss en el que se encuentra el DSSD_5
    int P_min=0; //Número de canal más bajo de los canales P
    int P_max=0; //Número de canal más alto de los canales P
    int N_min=0; //Número de canal más bajo de los canales N
    int N_max=0; //Número de canal más alto de los canales N
    
    if (numCalibDetector == 5){ //Parámetros para creal el histograma Raw del DSSD 5 ( el único que nos interesa de momento)
        nDSSD=4;
        n_buss=5;
        P_min=0;
        P_max=15;
        N_min=32; 
        N_max=47;
    }
    
    int nCh = 0; //Número de canal en el que se encuentra el detector.
    int nPAD = 0; // Número de histograma que corresponde al detector
    
    if (numCalibDetector == 5){//Parámetros para el PAD 5
        nCh = 5;
        nPAD = 4;
    }
    
    //Para evitar el error si no hay entries
    int printFrequency = (nEntries > 100) ? (nEntries / 100) : 1;

    //Bucle sobre las entradas
    for (int i = 0; i < int(nEntries); i++){
        //Cogemos el evento, super importante
        chain-> GetEvent(i);
        if (CalibDetector == 2 || CalibDetector == 0){
            for (int j = 0; j < vmmr8_counter; j++){// Bucle por si hubiera más VMMR8
                for (Int_t k = 0; k < mul_VMMR8[j]; k++){
                    
                    if (buss_VMMR8[j][k] != n_buss) continue; //Si no se encuentra en el buss al que está conectado el DSSD que nos interesa pasa al siguiente.
                    
                    int idx= detCh_VMMR8[j][k];
                    if (idx < P_min || idx > N_max) continue;
                    if (idx > P_max && idx < N_min) continue; //Nos quedamos sólo con los canales con strips de los DSSD conectadas
                    
                    if (idx > P_max){ //Evitamos irnos fuera del rango del array (max 31)
                        idx=idx-16; 
                    }
                    DSSD_Raw_hist_Total[nDSSD][idx]->Fill(values_VMMR8[j][k]);//Values son los canales electrónicos, detCh los físicos (detectores)
                }
            }
        }
        
        if (CalibDetector == 3 || CalibDetector == 0){
            for (int j = 0; j < mdpp32_counter; j++){// Bucle por si hubiera mas MDPP32
                for (Int_t k = 0; k < mul_MDPP32[j]; k++){ 
                    int idx= detCh_MDPP32[j][k];
                    if (idx != nCh) continue; //Nos quedamos sólo con eñ PAD que nos interesa
                    PAD_Raw_hist_Total[nPAD]->Fill(values_MDPP32[j][k]);//Values son los canales electrónicos, detCh los físicos (detectores)
                }
                if (j == 1) continue;
            }
        }
            
        if (i % printFrequency == 0 || i == nEntries - 1) {  // Imprimir cada 1% o en la última iteración
            double progress = (double(i) / nEntries) * 100.0;
            std::cout << " Creando histos Raw: " << progress << "% c.\r";
            std::cout.flush();  // Asegurar que la salida se actualice en la consola
        }
    }
    
    if (CalibDetector == 2){
        
        //Es muy difícil que TSpectrum sea capaz de identificar los picos que se solapan totalmente, por lo que vamos a utilizar para que los busque los 4 picos que se ven claros y a partir de ahí se hará un ajuste a dos gaussianas.
        
        std::vector<double> Calibration_Values_DSSD =Select_Calibration_Values(isotope);
        std::vector<double> Ranges_values_Left = Select_Calibration_Ranges("Alpha DSSD Left",CalibDetector);
        std::vector<double> Ranges_values_Right = Select_Calibration_Ranges("Alpha DSSD Right",CalibDetector);
        const int nPeaks_DSSD = 4;

        //Se crean los vectores en los que se van a guardar los datos
        std::vector<double> centroides_raw_DSSD(Calibration_Values_DSSD.size(),0.0);
        std::vector<double> err_centroides_DSSD(Calibration_Values_DSSD.size(),0.0);
        
        //Se crean las variables de calibración
        std::vector<double> a(32.0,0.0), err_a(32.0,0.0);
        std::vector<double> b(32.0,0.0), err_b(32.0,0.0);
        std::vector<int> entradas_dig(32), canales_fis(32);
                
        //Se crea un noise_cut para evitar el ruido de bajas energías, xmax corresponderá al canal máximo en el que se buscará
        int xmax_DSSD=2000;
        double sigma_tentativa=10;
        
        if (numCalibDetector != 0){ //Si sólo se ha pedido calibrar sólo uno de los detectores
            
            //Bucle que recorre todas las strips del DSSD para hacer la calibración de cada una de ellas
            for (int i = 0; i < 32; i++){
                
                //Se crea el histograma sobre el que se van a aplicar las funciones
                TH1F* histo_DSSD = DSSD_Raw_hist_Total[nDSSD][i];
                
                //Se crea la variable peaks_DSSD, que es donde se guardarán los picos encontrados por la función Buscar_Picos_DSSD
                std::vector<double> peaks_DSSD;
                Buscar_Picos_DetSil(histo_DSSD, noise_cut_dig[4][i], peaks_DSSD, nPeaks_DSSD, isotope, xmax_DSSD,sigma_tentativa);

                //Se limpia las variables de los centroides, para evitar errores en los bucles
                centroides_raw_DSSD.clear();
                err_centroides_DSSD.clear();                    
                double sigma_estimada;//Va a utilizar la media obtenida de la primera Gaussiana como sigma estimada para el resto

                
                //Se va a hacer el ajuste pico a pico
                for (int j = 0; j < peaks_DSSD.size(); j++){
                    double peakValue = peaks_DSSD[j];

                    double mean, sigma, area;
                    double err_mean;
                    
                    //Se hace el ajuste para obtener el centroide en canal
                    Ajuste_Gauss(histo_DSSD, peakValue, Ranges_values_Left[j], Ranges_values_Right[j], sigma_tentativa, area, mean, sigma, err_mean);//Ajuste Gaussiano
                        
                    centroides_raw_DSSD.push_back(mean);
                    err_centroides_DSSD.push_back(err_mean);            
                    sigma_estimada = sigma;
                    
                    double amplitud1, amplitud2, mean1, sigma1, mean2, sigma2, area1, area2, err_amplitud1, err_amplitud2, err_mean1, err_mean2, err_sigma1, err_sigma2, err_area1, err_area2;   
                }
                
                //Parámetros en los que se van a guardar los ajustes
                double a_lineal, b_lineal, err_a_lineal, err_b_lineal;
                    
                //Hacemos un ajuste lineal entre los centroides en canales encontrados con los valores de la energía tabulados. Con esto se obtienen los parámetros de calibración.
                Ajuste_Lineal(centroides_raw_DSSD, Calibration_Values_DSSD, err_centroides_DSSD, a_lineal, b_lineal, err_a_lineal, err_b_lineal);
                    
                //Guardamos en el vector las calibraciones obtenidas para cada strip
                a[i] = a_lineal;
                b[i] = b_lineal;
                err_a[i] = err_a_lineal;
                err_b[i] = err_b_lineal;
                
                //Guardamos a qué canal electrónico y físico pertenecen
                canales_fis[i] = i;
                if (i < 16){
                    entradas_dig[i]=i;
                }
                else {
                    entradas_dig[i]=i+16;
                }                
            double progress = (double(i) / 32) * 100.0;
            std::cout << " Realizando calibraciones del DSSD"<< nDSSD<<": " << progress << "% c.\r";
            std::cout.flush();  // Asegurar que la salida se actualice en la consola
            }
            //Extraer y guardar calibración
            std::string ArchCalib = "Calibracion/Calibration_DSSD_" + filename + ".txt";
            
            // Guardar calibración manual en txt
            ActualizarCalibracionDSSD(ArchCalib, nDSSD+1, n_buss, entradas_dig, canales_fis, a,b); 
        }
    }
    
    else if(CalibDetector == 3){
        
        //Es muy difícil que TSpectrum sea capaz de identificar los picos que se solapan totalmente, por lo que vamos a utilizar para que los busque los 4 picos que se ven claros y a partir de ahí se hará un ajuste a dos gaussianas.
        
        std::vector<double> Calibration_Values_PAD=Select_Calibration_Values(isotope);
        std::vector<double> Ranges_values = Select_Calibration_Ranges(isotope,CalibDetector);
        const int nPeaks_PAD = 4;

        //Se crean los vectores en los que se van a guardar los datos
        std::vector<double> centroides_raw_PAD(Calibration_Values_PAD.size(),0.0);
        std::vector<double> err_centroides_PAD(Calibration_Values_PAD.size(),0.0);
        
        double a_PAD = 0.0, err_a_PAD= 0.0;
        double b_PAD = 0.0, err_b_PAD = 0.0;
        int entrada_dig_PAD = nCh;
        int canal_fis_PAD = nPAD + 1;
        
        //Se crea un noise_cut para evitar el ruido de bajas energías
        double noise_cut_PAD = 5000;
        int xmax_PAD=20000;
        double sigma_tentativa = 100;
        
        if (numCalibDetector != 0){ //Si sólo se ha pedido calibrar sólo uno de los detectores
            
            //Se crea el histograma sobre el que se van a aplicar las funciones
            TH1F* histo_PAD = PAD_Raw_hist_Total[nPAD];
            
            //Se crea la variable peaks_PAD, que es donde se guardarán los picos encontrados por la función Buscar_Picos_DetSil
            std::vector<double> peaks_PAD;
            Buscar_Picos_DetSil(histo_PAD, noise_cut_PAD, peaks_PAD, nPeaks_PAD, isotope, xmax_PAD, sigma_tentativa);
            
            //Se limpia las variables de los centroides, para evitar errores en los bucles
            centroides_raw_PAD.clear();
            err_centroides_PAD.clear();
                
            //Se va a hacer el ajuste pico a pico
            for (int j = 0; j < peaks_PAD.size(); j++){
                double peakValue = peaks_PAD[j];

                double mean, sigma, area;
                double err_mean;
                    
                //Se hace el ajuste para obtener el centroide en canal
                Ajuste_Gauss(histo_PAD, peakValue, Ranges_values[j], Ranges_values[j],sigma_tentativa, area, mean, sigma, err_mean);//Ajuste Gaussiano
                    
                centroides_raw_PAD.push_back(mean);
                err_centroides_PAD.push_back(err_mean);            
            }
                                    
            //Hacemos un ajuste lineal entre los centroides en canales encontrados con los valores de la energía tabulados. Con esto se obtienen los parámetros de calibración.
            Ajuste_Lineal(centroides_raw_PAD, Calibration_Values_PAD, err_centroides_PAD, a_PAD, b_PAD, err_a_PAD, err_b_PAD);
        }
            
            
        //Extraer y guardar calibración
        std::string ArchCalib = "Calibracion/Calibration_PAD_" + filename + ".txt";
        
        // Guardar calibración manual en txt
        ActualizarCalibracionPAD(ArchCalib, entrada_dig_PAD, canal_fis_PAD, a_PAD,b_PAD); 
    }
    
    //Guardamos los histogramas creados
    gROOT->SetBatch(kTRUE);

    TDirectory* DSSDsDir = outputFile->mkdir("DSSDs");
    Save_DSSD_Histos(DSSDsDir, "Raw", DSSD_Raw_hist_Total);

    TDirectory* PADsDir = outputFile->mkdir("PADs");
    Save_PAD_Histos(PADsDir, "Raw", PAD_Raw_hist_Total);
}
