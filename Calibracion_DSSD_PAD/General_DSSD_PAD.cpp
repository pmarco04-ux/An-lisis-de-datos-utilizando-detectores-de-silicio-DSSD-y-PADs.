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
#include <iomanip>

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

#include <stdexcept>
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMatrixD.h"


//Implementación TSpectrum

void Buscar_Picos_DetSil(TH1F* histo, double noise_cut, std::vector<double> &peaks, const int nPeaks, std::string isotope, const int xmax, const double sigma_esp) {
    //Canal máximo hasta el que se va a buscar (xmax)
    gErrorIgnoreLevel = kFatal;
    
    // Se llama la función TSpectrum para que actúe sobre el número de picos de la fuente
    TSpectrum spectrum(nPeaks);

    // Se limita el rango de canales en los que se van a buscar los picos
    histo->GetXaxis()->SetRangeUser(noise_cut, xmax);
    
    //En nfound se guarda el número de picos encontrados. La función .Search tiene como parámetros el histograma en el que se va a buscar, la sigma esperada (20), alguna configuración opcional, y la altura mínima de los picos respecto al máximo del histograma (20%)
    int nfound = spectrum.Search(histo, sigma_esp, "", 0.2);
    
    //A continuación se sacan las posiciones de los picos y se ordenan de mayor a menor
    double *xpeaks = spectrum.GetPositionX();
    peaks.clear();
    for (int i = 0; i < nfound; i++) {
        peaks.push_back(xpeaks[i]);
    }
    std::sort(peaks.begin(), peaks.end());
    histo->GetXaxis()->SetRangeUser(0, xmax); // restaurar rango original para el dibujo
    gErrorIgnoreLevel = kInfo;
}

//Ajustes

//Ajuste lineal, importante para hacer la calibración

void Ajuste_Lineal(std::vector<double> x, std::vector<double> y, std::vector<double> sigma_y, double &a, double &b, double &err_a, double &err_b) {
    int n = (int)x.size();
    if (n != (int)y.size() || n != (int)sigma_y.size()) {
        std::cerr << "Ajuste_Lineal: vectores de distinto tamaño x: " << n << " Energias: "<< y.size() << " sigma " << sigma_y.size()<< std::endl;
        return;
    }

    double S=0, Sx=0, Sy=0, Sxx=0, Sxy=0;
    for (int i = 0; i < n; i++) {
        if (sigma_y[i] <= 0) {
            std::cerr << "Ajuste_Lineal: sigma <= 0 en punto " << i << ", se omite" << std::endl;
            continue;
        }
        double w = 1.0 / (sigma_y[i] * sigma_y[i]); // peso = 1/σ²
        S   += w;
        Sx  += w * x[i];
        Sy  += w * y[i];
        Sxx += w * x[i] * x[i];
        Sxy += w * x[i] * y[i];
    }

    double Delta = S * Sxx - Sx * Sx;
    if (Delta == 0) {
        std::cerr << "Ajuste_Lineal: Delta=0, ajuste imposible" << std::endl;
        return;
    }

    a     = (S   * Sxy - Sx * Sy)  / Delta;
    b     = (Sxx * Sy  - Sx * Sxy) / Delta;
    err_a = TMath::Sqrt(S   / Delta);
    err_b = TMath::Sqrt(Sxx / Delta);
}

void Ajuste_Gauss(TH1F* Histo, double media_tentativa, double rango_izq, double rango_der, double sigma_tentativa, Double_t &area, Double_t &mean, Double_t &sigma, Double_t &err_mean){
    
    //Se crea una función gauss a la que se va a ajustar dentro del rango dado
    
    TF1 *gauss = new TF1("gauss", "gaus", media_tentativa-rango_izq, media_tentativa+rango_der); //Defino la funcion de ajuste y el rango a analizar y donde
    gauss->SetLineColor(kRed); //Color del ajuste


    //Inicializo parámetros
    gauss->SetParameter(0, Histo->GetMaximum()); //Amplitud
    gauss->SetParameter(1, media_tentativa); //Media
    gauss->SetParameter(2, sigma_tentativa); //Sigma
    gauss->SetParLimits(2,0.01,1e9); //Límites de sigma


    Histo->Fit(gauss,"RQ+");// Ajusto mi histograma, el RQ+ hace que se vea laf
    
    //Cojo los parametros
    Double_t A = gauss->GetParameter(0);//Extrae la amplitud
    mean = gauss->GetParameter(1);//Extra la media
    err_mean = gauss-> GetParError(1);//Extra el error de la media
    sigma = gauss->GetParameter(2);//Extra la sigma
    double bin_width =  Histo->GetBinWidth(1);//Extrae la sigma
    area = A * TMath::Sqrt(2 * TMath::Pi()) * sigma/bin_width; //Fórmula del área gaussiana
    
    delete gauss; //Se borra para no crear fugas de memoria
}

//Introducción al ajuste a dos gaussianas. Será importante porque en la mayoría de los picos se solapan dos picos de forma que se pueden separar.
void Ajuste_Dos_Gauss(TH1F* Histo, double media_tentativa2, double rango_izq, double rango_der,double sigma_estimada, double intensidad_rel, double &amplitud1, double &amplitud2, double &mean1, double &mean2, double &sigma1, double &sigma2, double &area1, double &area2, double &err_amplitud1, double &err_amplitud2, double &err_mean1, double &err_mean2, double &err_sigma1, double &err_sigma2, double &err_area1, double &err_area2){

    //Para definir la segunda media tentativa utilizo la sigma estimada
    double media_tentativa1 = media_tentativa2-0.7*sigma_estimada;
    
    //Defino el rango entre el que se va a hacer el ajuste
    double xmin = media_tentativa1-rango_izq;
    double xmax = media_tentativa2+rango_der;
    
    //Creamos la función a la que se va a ajustar
    TF1 *gauss_gauss = new TF1("gauss_gauss", "gaus(0) + gaus(3)", xmin, xmax);
    
    //Definimos las alturas esperadas de cada uno de los picos. Para ello utilizamos el 
    double yPico2 = Histo->GetBinContent(Histo->FindBin(media_tentativa2));
    double yPico1 = yPico2*intensidad_rel;
    
    //Doy unos parámetros iniciales
    gauss_gauss->SetLineColor(kBlack);
    gauss_gauss->SetParameter(0, yPico1);   // amplitud gaussiana
    gauss_gauss->SetParameter(1, media_tentativa1); // media
    gauss_gauss->SetParameter(2, sigma_estimada); //sigma inicial (ajusta según tu resolución)                 
    gauss_gauss->SetParameter(3, yPico2); // amplitud gaussiana
    gauss_gauss->SetParameter(4, media_tentativa2); // media
    gauss_gauss->SetParameter(5, sigma_estimada); //sigma inicial (ajusta según tu resolución)
    
    Histo->Fit(gauss_gauss, "RQ+");
    amplitud1  = gauss_gauss->GetParameter(0);
    mean1        = gauss_gauss->GetParameter(1);
    sigma1       =abs(gauss_gauss->GetParameter(2));
    amplitud2  = gauss_gauss->GetParameter(3);
    mean2        = gauss_gauss->GetParameter(4);
    sigma2       =abs(gauss_gauss->GetParameter(5));
    
    err_amplitud1 = gauss_gauss->GetParError(0);
    err_mean1     = gauss_gauss->GetParError(1);
    err_sigma1    = gauss_gauss->GetParError(2);
    err_amplitud2 = gauss_gauss->GetParError(3);
    err_mean2     = gauss_gauss->GetParError(4);
    err_sigma2    = gauss_gauss->GetParError(5);
    
    double bin_width =  Histo->GetBinWidth(1);
    area1 = amplitud1 * TMath::Sqrt(2 * TMath::Pi()) * sigma1/bin_width; //Habrá que dividirla entre el width de los canales
    area2 = amplitud2 * TMath::Sqrt(2 * TMath::Pi()) * sigma2/bin_width; //Habrá que dividirla entre el width de los canales
    err_area1 = (TMath::Sqrt(2 * TMath::Pi()) *TMath::Sqrt(TMath::Power(sigma1 * err_amplitud1, 2) +TMath::Power(amplitud1 * err_sigma1, 2)))/bin_width;
    err_area2 = (TMath::Sqrt(2 * TMath::Pi()) *TMath::Sqrt(TMath::Power(sigma2 * err_amplitud2, 2) +TMath::Power(amplitud2 * err_sigma2, 2)))/bin_width;

    //Vamos a graficar ambas gaussianas por separado
    TF1 *g1 = new TF1("g1", "gaus", xmin, xmax);
    TF1 *g2 = new TF1("g2", "gaus", xmin, xmax);
    g1->SetParameters(amplitud1, mean1, sigma1);
    g2->SetParameters(amplitud2, mean2, sigma2);
    g1->SetLineColor(kBlue);
    g1->SetLineStyle(2);
    g2->SetLineColor(kRed);
    g2->SetLineStyle(2);

    // Las dejamos enganchadas al histograma para que se vean también si lo abres más tarde
    Histo->GetListOfFunctions()->Add(g1);
    Histo->GetListOfFunctions()->Add(g2);

    // Dibujo en directo
    static TCanvas* c_ajuste = nullptr;
    if (!c_ajuste) c_ajuste = new TCanvas("c_ajuste_dos_gauss", "Ajuste doble gaussiana", 800, 600);
    c_ajuste->cd();
    Histo->Draw();       // dibuja el histo + gauss_gauss (ya enganchado por el Fit con "+")
    g1->Draw("same");
    g2->Draw("same");
    c_ajuste->Update();
    gSystem->ProcessEvents();
}

void Save_Histo(TObject* Histo, TDirectory* Directory_name) {
    Directory_name->cd();
    Histo->Write();
    delete Histo;
}


std::vector <double> Select_Calibration_Values(const std::string &isotope){
    
    if(isotope == "Alpha PAD"){
        return {3182.69, 5148.89, 5480.023, 5794.88};//Fuente combinada de Gd, Pu, Am, Cm
    }
    else if(isotope == "Alpha DSSD"){
        return {3182.69, 5148.89, 5480.023, 5794.88};//Fuente combinada
    }
    
    else{
        throw std::runtime_error("Isotope not found: " + isotope);
    }
}

std::vector <double> Select_Calibration_Ranges(const std::string &isotope,const int &CalibDetector){

    if (isotope == "Alpha PAD"){
        return {500,500,500,500};
    }
    
    else if (isotope == "Alpha DSSD Left"){
        return {9,10,10,10};
    }
    
        else if (isotope == "Alpha DSSD Right"){
        return {14,14,14,14};
    }
    
    else{
        throw std::runtime_error("Isotope not found: " + isotope);
    }
}

std::vector <double> Select_Fit_Ranges(const std::string &isotope){

    if(isotope == "Alpha PAD"){
        return {30,30,30,30};
    }
    
    else if(isotope == "Alpha DSSD"){
        return {30,30,30,30};
    }
    
    else{
        throw std::runtime_error("Isotope not found: " + isotope);
    }
}

void Draw_Calib_lines(TH1F* Histo, const double Calib_value, double height=800,int line_color= kRed){
    double x_Points[2] = {Calib_value,Calib_value};
    double y_Points[2] = {0.0,height};
    TPolyLine* Calib_line = new TPolyLine(2,x_Points,y_Points);
    Calib_line ->SetLineColor(line_color);
    Calib_line ->SetLineWidth(2);
    Calib_line ->SetLineStyle(2);
    Histo->GetListOfFunctions()->Add(Calib_line); 
}

void Draw_Diagonal_2D(TH2F* Histo, const double min, const double max, int line_color = kRed){
    double x_Points[2] = {min,max};
    double y_Points[2] = {min,max};
    TPolyLine* Diagonal_line = new TPolyLine(2,x_Points,y_Points);
    Diagonal_line->SetLineColor(line_color);
    Diagonal_line->SetLineWidth(2);
    Diagonal_line->SetLineStyle(2);
    Histo->GetListOfFunctions()->Add(Diagonal_line);
}

void Draw_Check_Calib_lines(TH2F* Histo, const double Calib_value){
    //Toma como límites los límites del histograma
    double x_min = Histo->GetXaxis()->GetXmin();
    double x_max = Histo->GetXaxis()->GetXmax();
    
    double x_Points[2] = {x_min,x_max};
    double y_Points[2] = {Calib_value,Calib_value};
    TPolyLine* Calib_line = new TPolyLine(2,x_Points,y_Points);
    Calib_line ->SetLineColor(kRed);
    Calib_line ->SetLineWidth(2);
    Calib_line ->SetLineStyle(2);
    Histo->GetListOfFunctions()->Add(Calib_line); 
}


void Save_BKG(TFile* outputFile, 
    TH1F* GAGGS_BKG_hist_Long[8],
    TH1F* GAGGS_BKG_hist_Short[8],
    TH1F* GAGGS_BKG_Total_Long,
    TH1F* GAGGS_BKG_Total_Short
    ){
    gROOT->SetBatch(kTRUE);
    TDirectory* BKG_Dir =  outputFile->mkdir("BKG");
    TDirectory* GAGG_Calib_Long_Dir = BKG_Dir->mkdir("Calibrated Integration Long");
    TDirectory* GAGG_Calib_Short_Dir = BKG_Dir->mkdir("Calibrated Integration Short");
    TDirectory* Total_Dir = BKG_Dir->mkdir("Total");

    for(int i = 0; i < 8; i ++){
        Save_Histo(GAGGS_BKG_hist_Long[i],GAGG_Calib_Long_Dir);
        Save_Histo(GAGGS_BKG_hist_Short[i],GAGG_Calib_Short_Dir);
    }
    Save_Histo(GAGGS_BKG_Total_Long,Total_Dir);
    Save_Histo(GAGGS_BKG_Total_Short,Total_Dir);
    std::cout << " BKG guardado " << "\n";
}

//Esta función sirve para reescribir sólo una parte de los archivos de calibración, para que al calibrar un sólo detector no se pierda el resto. Para ello, primero vamos a abrir el documento para leerlo y analizarlo antes de reescribir

struct CalibParam_DSSD {
    int modulo;
    int buss;
    int entrada_digital;
    int canal_fisico;
    double a;
    double b;
};

void ActualizarCalibracionDSSD(
    const std::string& filename, 
    int mod_target, 
    int buss_target, 
    const std::vector<int>& entradas_dig, 
    const std::vector<int>& canales_fis, 
    const std::vector<double>& nuevas_a, 
    const std::vector<double>& nuevas_b) 
{
    std::vector<CalibParam_DSSD> parametros_guardados;

    // 1. Leer el archivo si existe y conservar lo que NO sea del MÓDULO y BUSS actual
    std::ifstream fileIn(filename);
    if (fileIn.is_open()) {
        std::string line;
        while (std::getline(fileIn, line)) {
            if (line.empty() || line[0] == '#') continue; // Ignorar comentarios

            std::stringstream ss(line);
            CalibParam_DSSD p;
            if (ss >> p.modulo >> p.buss >> p.entrada_digital >> p.canal_fisico >> p.a >> p.b) {
                // Conservar solo si pertenece a OTRO módulo o a OTRO buss
                if (p.modulo != mod_target || p.buss != buss_target) {
                    parametros_guardados.push_back(p);
                }
            }
        }
        fileIn.close();
    }

    // 2. Añadir las nuevas calibraciones para este MÓDULO y BUSS
    for (size_t i = 0; i < nuevas_a.size(); ++i) {
        CalibParam_DSSD p;
        p.modulo = mod_target;
        p.buss = buss_target;
        p.entrada_digital = entradas_dig[i];
        p.canal_fisico = canales_fis[i];
        p.a = nuevas_a[i];
        p.b = nuevas_b[i];
        parametros_guardados.push_back(p);
    }

    // 3. Sobrescribir el archivo conservando la cabecera y el formato de 6 columnas
    std::ofstream fileOut(filename);
    if (!fileOut.is_open()) {
        std::cerr << "Error al abrir el archivo de salida: " << filename << std::endl;
        return;
    }

    fileOut << "# DSSD   BUSS   ENTRADA_DIGITAL   CANAL_FISICO   Slope_a       Intercept_b\n";
    for (const auto& p : parametros_guardados) {
        fileOut << std::fixed << std::setprecision(8)
                << std::setw(8)  << p.modulo
                << std::setw(8)  << p.buss
                << std::setw(18) << p.entrada_digital
                << std::setw(15) << p.canal_fisico
                << std::setw(16) << p.a
                << std::setw(16) << p.b << "\n";
    }
    fileOut.close();
    std::cout << "Calibración actualizada correctamente en: " << filename << std::endl;
}

struct CalibParam_PAD {
    int entrada_digital;
    int canal_fisico;
    double a;
    double b;
};


void ActualizarCalibracionPAD(const std::string& filename, int entrada_dig, int canal_fis, double nueva_a, double nueva_b) 
{
    std::vector<CalibParam_PAD> parametros_guardados;

    // 1. Leer el archivo si existe y conservar los PADs distintos al actual
    std::ifstream fileIn(filename);
    if (fileIn.is_open()) {
        std::string line;
        while (std::getline(fileIn, line)) {
            if (line.empty() || line[0] == '#') continue; // Ignorar comentarios

            std::stringstream ss(line);
            CalibParam_PAD p;
            if (ss >> p.canal_fisico >> p.entrada_digital >> p.a >> p.b) {
                // Conservar solo si pertenece a OTRO canal físico o digital
                if (p.canal_fisico != canal_fis || p.entrada_digital != entrada_dig) {
                    parametros_guardados.push_back(p);
                }
            }
        }
        fileIn.close();
    }

    // 2. Añadir la nueva calibración para este PAD (sin bucle for)
    CalibParam_PAD p;
    p.canal_fisico = canal_fis;
    p.entrada_digital = entrada_dig;
    p.a = nueva_a;
    p.b = nueva_b;
    parametros_guardados.push_back(p);

    // 3. Sobrescribir el archivo
    std::ofstream fileOut(filename);
    if (!fileOut.is_open()) {
        std::cerr << "Error al abrir el archivo de salida: " << filename << std::endl;
        return;
    }

    fileOut << "# PAD   ENTRADA_DIGITAL   Slope_a       Intercept_b\n";
    for (const auto& param : parametros_guardados) {
        fileOut << std::fixed << std::setprecision(8)
                << std::setw(15) << param.canal_fisico
                << std::setw(18) << param.entrada_digital
                << std::setw(16) << param.a
                << std::setw(16) << param.b << "\n";
    }
    fileOut.close();
    std::cout << "Calibración actualizada correctamente en: " << filename << std::endl;
}


//Cargar los parámetros de calibración del MDPP32
void LoadMDPP32CalibParameters(
    const std::string& inputCalibMDPP32,
    double (&a_MDPP32)[5],
    double (&b_MDPP32)[5])
{
    std::ifstream inputFile(inputCalibMDPP32);
    if (!inputFile.is_open()) {
        std::cerr << "Error al abrir el archivo de calibracion PAD: " << inputCalibMDPP32 << std::endl;
        return;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        if (line.empty() || line[0] == '#') continue; // Saltar comentarios y líneas vacías

        std::istringstream iss(line);
        int det = -1;
        int ent_dig = -1;
        double a_val = 0.0, b_val = 0.0;

        // Formato esperado según ActualizarCalibracionPAD:
        // canal_fisico (det) | entrada_digital | Slope_a | Intercept_b
        if (iss >> det >> ent_dig >> a_val >> b_val) {
            // Convertimos de base 1 (física) a base 0 (índice del array C++)
            int index = det - 1; 
            if (index >= 0 && index < 5) {
                a_MDPP32[index] = a_val;
                b_MDPP32[index] = b_val;
            }
        }
    }
    inputFile.close();
    std::cout << "Parámetros MDPP32 (PAD) cargados exitosamente." << std::endl;
}


//Cargar los parámetros de Calibración de los VMMR8
void LoadVMMR8CalibParameters(
    const std::string& inputCalibVMMR8,
    double (&a_VMMR8)[5][32],
    double (&b_VMMR8)[5][32])
{
    std::ifstream inputFile(inputCalibVMMR8);
    if (!inputFile.is_open()) {
        std::cerr << "Error al abrir el archivo de calibracion DSSD: " << inputCalibVMMR8 << std::endl;
        return;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        if (line.empty() || line[0] == '#') continue; // Saltar comentarios y líneas vacías

        std::istringstream iss(line);
        int mod = -1, buss = -1, ent_dig = -1, strip = -1;
        double a_val = 0.0, b_val = 0.0;

        // Formato esperado según ActualizarCalibracionDSSD:
        // DSSD (mod) | BUSS | ENTRADA_DIGITAL | CANAL_FISICO (strip) | Slope_a | Intercept_b
        if (iss >> mod >> buss >> ent_dig >> strip >> a_val >> b_val) {
            int modIdx = mod - 1;     // Detector 1..5 -> Índice 0..4
            int stripIdx = strip ; // Strip 1..32  -> Índice 0..31

            if (modIdx >= 0 && modIdx < 5 && stripIdx >= 0 && stripIdx < 32) {
                a_VMMR8[modIdx][stripIdx] = a_val;
                b_VMMR8[modIdx][stripIdx] = b_val;
            }
        }
    }
    inputFile.close();
    std::cout << "Parámetros VMMR8 (DSSD) cargados exitosamente." << std::endl;
}

//Cargar los parámetros de Mapping de los DSSD
void LoadDSSDMappingParameters(
    const std::string& inputMappingDSSD,
    int (&Ch_fis)[5][32])
{
    std::ifstream inputFile(inputMappingDSSD);
    if (!inputFile.is_open()) {
        std::cerr << "Error al abrir el archivo de mapeado DSSD: " << inputMappingDSSD << std::endl;
        return;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        if (line.empty() || line[0] == '#') continue; // Saltar comentarios y líneas vacías

        std::istringstream iss(line);
        int dssd = -1, ent_dig = -1, strip = -1;
        std::string tipo = "";

        // Formato esperado según Mapping_DSSD_IS690B:
        // DSSD | TIPO | ENTRADA_DIGITAL | CANAL_FISICO (strip) 
        if (iss >> dssd >> tipo >> ent_dig >> strip) {
            int ndssd = dssd - 1;     // Detector 1..5 -> Índice 0..4
            int Ch_dig = ent_dig ; // Canal digital

            if (ndssd >= 0 && ndssd < 5 && Ch_dig >= 0 && Ch_dig <= 47) {
                if (Ch_dig <= 15){
                    Ch_fis[ndssd][Ch_dig] = strip;
                }
                else if (Ch_dig >=32){
                    Ch_fis[ndssd][Ch_dig-16] = strip;
                }
            }
        }
    }
    inputFile.close();
    std::cout << "Mapping (DSSD) cargado exitosamente." << std::endl;
}

void LoadDSSDNoise_CutRawParameters(
    const std::string& inputNoise_CutRaw,
    int (&noise_cut_dig)[5][32], int (&noise_cut_fis)[5][32])
{
    std::ifstream inputFile(inputNoise_CutRaw);
    if (!inputFile.is_open()) {
        std::cerr << "Error al abrir el archivo de noise_cut DSSD: " << inputNoise_CutRaw<< std::endl;
        return;
    }

    std::string line;
    while (std::getline(inputFile, line)) {
        if (line.empty() || line[0] == '#') continue; // Saltar comentarios y líneas vacías

        std::istringstream iss(line);
        int dssd = -1, ch = -1, NC_fis = -1, NC_dig;
        std::string tipo = "";

        // Formato esperado según Mapping_DSSD_IS690B:
        // DSSD | CANAL | NOISE CUT DIGITAL | NOISE CUT FISICO 
        if (iss >> dssd >> ch >> NC_dig >> NC_fis) {

            if (dssd > 0 && dssd <= 5 && ch > 0 && ch <= 32) {
                noise_cut_fis[dssd-1][ch-1] = NC_fis;
                noise_cut_dig[dssd-1][ch-1] = NC_dig;
            }
        }
    }
    inputFile.close();
    std::cout << "Noise_cut (DSSD) cargado exitosamente." << std::endl;
}

    // Guarda los 5x32 histogramas de un DSSD (Raw o Calib) en outputFile,
    // organizados como DSSDs/<tipo>/DSSD_N/
    void Save_DSSD_Histos(TDirectory* DSSDsDir, const char* tipo, TH1F* hist[5][32]) {
        
        TDirectory* tipoDir = DSSDsDir->mkdir(tipo); 
        
        for (int det = 0; det < 5; det++) {
            TString detDirName = Form("DSSD_%d", det + 1);
            TDirectory* CurrentDetDir = tipoDir->mkdir(detDirName);
            CurrentDetDir->cd();
            
            for (int strip = 0; strip < 32; strip++) {
                Save_Histo(hist[det][strip], CurrentDetDir);
            }
        }
    }
    
    void Save_DSSD_Histos_2D(TDirectory* DSSDsDir, const char* tipo, TH2F* hist[5][32]) {
        TDirectory* tipoDir = DSSDsDir->mkdir(tipo);
        for (int det = 0; det < 5; det++) {
            TDirectory* CurrentDetDir = tipoDir->mkdir(Form("DSSD_%d", det + 1));
            for (int strip = 0; strip < 32; strip++) {
                Save_Histo(hist[det][strip], CurrentDetDir);
            }
        }
    }   

    // Guarda los 5 histogramas de PAD (Raw o Calib) en outputFile,
    // organizados como PADs/<tipo>/
    void Save_PAD_Histos(TDirectory* PADsDir, const char* tipo, TH1F* hist[5]) {
        
        TDirectory* tipoDir = PADsDir->mkdir(tipo); 
        
        for (int i = 0; i < 5; i++) {
            Save_Histo(hist[i], tipoDir);
        }
    }
    

int General_DSSD_PAD(){ 
    return 0;
}

