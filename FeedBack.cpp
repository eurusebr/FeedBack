// Here is a program to caclulate the feedback model
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>
#include <sys/stat.h>
#include <boost/numeric/odeint.hpp>
#include <chrono>
#include <boost/random.hpp>

using namespace std;
using namespace boost::numeric::odeint;
using namespace std::chrono;

// definition of imaginary unit
const complex<double> I(0.0, 1.0);

// defining functions
void GaussianArray(vector<complex<double>> &complexArray, double mean, double sigma);
double c_array_norm(vector<complex<double>> &a, int N);
vector<complex<double>> HDNLS(vector<complex<double>> &a, double alpha, double dt, int N);
void noise(vector<complex<double>> &a, vector<complex<double>> &b, boost::random::mt19937 &rng, boost::normal_distribution<> &m_dist, double gamma, double beta, double dt, int N);
void normalize(vector<complex<double>> &complexArray, int N);
void localizedArray(vector<complex<double>> &complexArray);
void normarray(vector<vector<double>> &a, vector<complex<double>> b, int row, double ens);
//void write(vector<double> &A, string &file);
void write(const vector<double> &data, const string &file);
void writeline(const vector<double> &data, const string &file);
vector<complex<double>> Enoughtpart(vector<complex<double>> &a, double alpha, double Enought, double newnorm, double xi, double lambda, double dt, int N);
complex<double> H(vector<complex<double>> &a, double alpha, int N);
vector<complex<double>> EnoughtpartDelta(vector<complex<double>> &a, double alpha, double Enought, double newnorm, double xi, double lambda, double dt, int N, double delta, int bond, int place);
void writecomp(const std::string& file_path, vector<complex<double>>& complex_vector);
void readvec(const std::string& filePath, vector<complex<double>>& complexVector, size_t N);
void NormalizedGaussian(std::vector<std::complex<double>>& result, double peak, double standardDeviation);
void writecompadd(const std::string& file_path, vector<complex<double>>& complex_vector);
void saveComplexVector(const std::string& file_path, const std::vector<std::complex<double>>& complex_vector);
std::vector<std::complex<double>> readComplexVector(const std::string& file_path);
int main(int argc, char **argv)
{

#include <iostream>
#include <vector>
#include <sstream>

    if (argc < 14) {
        cerr << "Usage: " << argv[0] << " <model> <other_arguments>" << endl;
        return 1;
    }

    string model = argv[1];
    int N = stod(argv[2]);
    double sigma = stod(argv[3]);
    double dt = stod(argv[4]);
    double alpha = stod(argv[5]);
    double gamma = stod(argv[6]);
    double beta = stod(argv[7]);
    double lambda = stod(argv[8]);
    //double Enought = stod(argv[9]);
    double xi = stod(argv[9]);
    int steps = stod(argv[10]);
    double ens = stod(argv[11]);
    int samples = stod(argv[12]);
    double delta = stod(argv[13]); 
    int bond = stod(argv[14]);

    auto start = steady_clock::now();
    
    string file, file1, file2, file3, file4, file5;
    //vector<vector<double>> amplitude(steps, vector<double>(N, 0));
    vector<double> pzero, energy, ave;
    vector<double> subamp;
    vector<double> subbamp(N);
    vector<double> newnorm;
    vector<double> amp;

    boost::random::mt19937 rng;
    boost::normal_distribution<> m_dist(0.0, sqrt(0.5 * dt));
    vector<complex<double>> complexArray(N);
    vector<complex<double>> test(N);
    vector<complex<double>> noiseArray(N);
    vector<complex<double>> AddingArray(N);

    // cout<<"complex:"<<endl;
    // for (size_t s = 0; s < N; s++)
    // {
    //     std::cout << complexArray[s]<<" ";
    // }
    // std::cout << endl;

    if (model == "DNLS")
    {
        ofstream myfile("./DNLS/" + string(argv[5]) + "/log.txt",std::ofstream::trunc); 

        for (int l = 0; l < 15; l++)
        {
            myfile << argv[l] << endl;
        }
        myfile.close();
        // Create a new directory
        string folder = "./DNLS/" + string(argv[5]) + "/";
        struct stat info;

        if (stat(folder.c_str(), &info) != 0)
        {
            mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            file = folder + "/" + argv[5] + ".bin";
            file1 = folder + "/energy.bin";
            file2 = folder + "/pzero.bin";
            file3 = folder + "/norm.bin";
            file4 = folder + "/startvec.txt";
        }

        if (stat(folder.c_str(), &info) == 0)
        {
            file = folder + "/" + argv[5] + ".bin";
            ofstream outFile(file, ios::binary | ios::trunc);
            file1 = folder + "/energy.bin";
            file2 = folder + "/pzero.bin";
            file3 = folder + "/norm.bin";
            file4 = folder + "/startvec.txt";
            ofstream outFile1(file1, ios::binary | ios::trunc);
            ofstream outFile2(file2, ios::binary | ios::trunc);
            ofstream outFile3(file3, ios::binary | ios::trunc);
        }

        // GaussianArray(complexArray, 0.0, sigma);
        //localizedArray(complexArray);
        //normalize(complexArray, N);
        complexArray = readComplexVector(file4);
        std::cout << H(complexArray, alpha, N).real() << endl;
        newnorm.push_back(c_array_norm(complexArray, N));
        energy.push_back(H(complexArray, alpha, N).real());
        pzero.push_back(norm(complexArray[int(N / 2)]));
        //energy.push_back(H(complexArray, alpha, N).real());
        //pzero.push_back(amplitude[0][int(N / 2)]);
        //normarray(amplitude, complexArray, 0, ens);
        //std::memcpy(&subamp[0], &amplitude[0][int(N / 2) - 20], 40 * sizeof(double));
        //write(subamp, file);

        for (size_t j = 0; j < steps-1; j++)
        {
            int k = 0;
            complexArray = HDNLS(complexArray, alpha, dt, N);
            normalize(complexArray, N);
            //normarray(amplitude, complexArray, j+1, ens);
//            write(subamp, file);
            if (j % samples == 0)
            {
                k++;
                if (k% 1000 == 0)
                {
                    std::cout<<k<<endl;
                }
                energy.push_back(H(complexArray, alpha, N).real());
                pzero.push_back(norm(complexArray[int(N / 2)]));
                newnorm.push_back(c_array_norm(complexArray, N));
                //pzero.push_back(amplitude[j][int(N / 2)]);
                //std::memcpy(&subamp[0], &amplitude[j][int(N / 2) - 20], 40 * sizeof(double));
                //write(subamp, file);
            }
        }
        writeline(energy, file1);
        writeline(newnorm, file3);
        writeline(pzero,file2);
        //write(energy, file1);

    }

    if (model == "ConDNLS")
    {
        // Create a new directory
        string folder = "./ConDNLS(" + string(argv[5]) + ")";
        struct stat info;
        if (stat(folder.c_str(), &info) != 0)
        {
            mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            string secfolder = folder + "/xi(" + argv[9] + ")";
            mkdir(secfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            file = secfolder + "/lambda(" + argv[8] + ").bin";
            file1 = secfolder + "/energy.bin";
            file2 = secfolder + "/pzero.bin";
            file3 = secfolder + "/norm.bin";
            file4 = secfolder + "/startvec.bin";
        }
        else if (info.st_mode & S_IFDIR)
        {
            string secfolder = folder + "/xi(" + argv[9] + ")";
            if (stat(secfolder.c_str(), &info) != 0)
            {
                mkdir(secfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                file = secfolder + "/lambda(" + argv[8] + ").bin";
                file1 = secfolder + "/energy.bin";
                file2 = secfolder + "/pzero.bin";
                file3 = secfolder + "/norm.bin";
                file4 = secfolder + "/startvec.bin";
                ofstream outFile(file, ios::binary | ios::trunc);
                ofstream outFile1(file1, ios::binary | ios::trunc);
                ofstream outFile2(file2, ios::binary | ios::trunc);
                ofstream outFile3(file3, ios::binary | ios::trunc);
            }
            else
            {
                file = secfolder + "/lambda(" + argv[8] + ").bin";
                file1 = secfolder + "/energy.bin";
                file2 = secfolder + "/pzero.bin";
                file3 = secfolder + "/norm.bin";
                file4 = secfolder + "/startvec.bin";
                ofstream outFile(file, ios::binary | ios::trunc);
                ofstream outFile1(file1, ios::binary | ios::trunc);
                ofstream outFile2(file2, ios::binary | ios::trunc);
                ofstream outFile3(file3, ios::binary | ios::trunc);
            }
        }
        ofstream myfile("./ConDNLS(" + string(argv[5]) + ")/xi(" + string(argv[9]) + ")" + "/log.txt",std::ofstream::trunc); 
        for (int l = 0; l < 15; l++)
        {
            myfile << argv[l] << endl;
        }
        myfile.close();
        // GaussianArray(complexArray, 0.0, sigma);
        localizedArray(complexArray);
        normalize(complexArray, N);
        //complexArray = readComplexVector(file4);
        std::cout << "Energy = " << H(complexArray, alpha, N).real() << endl;
        newnorm.push_back(c_array_norm(complexArray, N));
        double normm = newnorm[0];
        energy.push_back(H(complexArray, alpha, N).real());
        pzero.push_back(norm(complexArray[int(N / 2)]));
        //energy.push_back(H(complexArray, alpha, N).real());
        //pzero.push_back(amplitude[0][int(N / 2)]);
        //normarray(amplitude, complexArray, 0, ens);
        for (int i = 480; i <= 520; ++i) 
        {
            subamp.push_back(std::abs(complexArray[i]));
        }
        subamp.clear();
        //std::memcpy(&subamp[0], &amplitude[0][int(N / 2) - 20], 40 * sizeof(double));
        //write(subamp, file);
        int k = 0;
        for (size_t j = 0; j < steps-1; j++)
        {
            //std::cout << j << endl;
            complexArray = Enoughtpart(complexArray, alpha, energy[0], normm, xi, lambda, dt, N);
            normm = c_array_norm(complexArray, N);
            //normalize(complexArray, N);
            //normarray(amplitude, complexArray, j+1, ens);
//            write(subamp, file);
            if (j % samples == 0)
            {   
                //std::cout << k <<endl;
                k++;
                if (k% 100 == 0)
                {
                    std::cout<<k<<endl;
                }
                energy.push_back(H(complexArray, alpha, N).real());
                pzero.push_back(norm(complexArray[int(N / 2)]));
                newnorm.push_back(normm);
                //pzero.push_back(amplitude[j][int(N / 2)]);
                for (int i = 480; i <= 520; ++i) 
                {
                    subamp.push_back(std::abs(complexArray[i]));
                }
                writeline(subamp, file);
                subamp.clear();
            }
        }
        writeline(energy, file1);
        writeline(newnorm, file3);
        writeline(pzero,file2);
        //write(energy, file1);

    }

    if (model == "FeedBack")
    {
        // Create a new directory
        string folder = "./FeedBack_beta(" + string(argv[7]) + ")";
        struct stat info;
        if (stat(folder.c_str(), &info) != 0)
        {
            mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            string secfolder = folder + "/alpha(" + argv[5] + ")";
            mkdir(secfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            file = secfolder + "/gamma(" + argv[6] + ").bin";
        }
        else if (info.st_mode & S_IFDIR)
        {
            string secfolder = folder + "/alpha(" + argv[5] + ")";
            if (stat(secfolder.c_str(), &info) != 0)
            {
                mkdir(secfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                file = secfolder + "/gamma(" + argv[6] + ").bin";
                ofstream outFile(file, ios::binary | ios::trunc);
            }
            else
            {
                file = secfolder + "/gamma(" + argv[6] + ").bin";
                ofstream outFile(file, ios::binary | ios::trunc);
            }
        }

        for (size_t s = 0; s < ens; s++)
        {
            std::cout << s << endl;
            localizedArray(complexArray);
            normalize(complexArray, N);
            //normarray(amplitude, complexArray, 0, ens);
            for (size_t i = 0; i < steps - 1; i++)
            {
                noise(complexArray, AddingArray, rng, m_dist, gamma, beta, dt, N);
                complexArray = HDNLS(complexArray, alpha, dt, N);
                for (size_t j = 0; j < N; j++)
                {
                    complexArray[j] = complexArray[j] + AddingArray[j];
                }
                normalize(complexArray, N);
                //normarray(amplitude, complexArray, i + 1, ens);
            }
        }

        std::cout << steps << endl;
        for (int j = 0; j < steps; j++)
        {
            // std::memcpy(&subamp[0], &amplitude[j][int(N/2) - 20], 40 * sizeof(double));
            // pzero.push_back(amplitude[j][int(N/2)]);
            // write(subamp, file);
            if (j % samples == 0)
            {
                //pzero.push_back(amplitude[j][int(N / 2)]);
            }
        }
        write(pzero, file);

        // write(pzero, file);
        // std::ofstream outFile(file, std::ios::binary | std::ios::app);
        // outFile.write((char*)pzero.data(), pzero.size() * sizeof(double));
        // outFile.close();
    }

    if (model == "Damping")
    {
        // Create a new directory
        string folder = "./Damping_alpha(" + string(argv[5]) + ")";
        struct stat info;
        if (stat(folder.c_str(), &info) != 0)
        {
            mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            string secfolder = folder + "/xi(" + argv[9] + ")";
            mkdir(secfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            file = secfolder + "/lambda(" + argv[8] + ").bin";
            file1 = secfolder + "/energy.bin";
            file2 = secfolder + "/pzero.bin";
            file3 = secfolder + "/norm.bin";
            file4 = secfolder + "/startvec.bin";
        }
        else if (info.st_mode & S_IFDIR)
        {
            string secfolder = folder + "/xi(" + argv[9] + ")";
            if (stat(secfolder.c_str(), &info) != 0)
            {
                mkdir(secfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                file = secfolder + "/lambda(" + argv[8] + ").bin";
                file1 = secfolder + "/energy.bin";
                file2 = secfolder + "/pzero.bin";
                file3 = secfolder + "/norm.bin";
                file4 = secfolder + "/startvec.bin";
                ofstream outFile(file, ios::binary | ios::trunc);
                ofstream outFile1(file1, ios::binary | ios::trunc);
                ofstream outFile2(file2, ios::binary | ios::trunc);
                ofstream outFile3(file3, ios::binary | ios::trunc);
            }
            else
            {
                file = secfolder + "/lambda(" + argv[8] + ").bin";
                file1 = secfolder + "/energy.bin";
                file2 = secfolder + "/pzero.bin";
                file3 = secfolder + "/norm.bin";
                file4 = secfolder + "/startvec.bin";
                ofstream outFile(file, ios::binary | ios::trunc);
                ofstream outFile1(file1, ios::binary | ios::trunc);
                ofstream outFile2(file2, ios::binary | ios::trunc);
                ofstream outFile3(file3, ios::binary | ios::trunc);
            }
        }


        for (size_t s = 0; s < ens; s++)
        {
            int k = 0;
            std::cout << s << endl;
            localizedArray(complexArray);
            //readvec(file4, complexArray, N);
            //NormalizedGaussian(complexArray, 0.828665, 0.8);
            newnorm.push_back(c_array_norm(complexArray, N));
            double normm = newnorm[0];
            //cout << c_array_norm(complexArray, N) << endl;
            //normalize(complexArray, N);
            energy.push_back(H(complexArray, alpha, N).real());
            //normarray(amplitude, complexArray, 0, ens);
            //std::memcpy(&subamp[0], &amplitude[0][int(N/2) - 20], 40 * sizeof(double));
            //write(subamp, file);
            for (size_t i = 0; i < steps - 1; i++)
            {
                //complexArray = Enoughtpart(complexArray, alpha, energy[0], newnorm, xi, lambda, dt, N);
                complexArray = EnoughtpartDelta(complexArray, alpha, energy[0], normm, xi, lambda, dt, N, delta, bond, i);
                normm = c_array_norm(complexArray, N);
                //normalize(complexArray, N);
                
                //normarray(amplitude, complexArray, i + 1, ens);
//                if (i % samples == 0)
//                {
//                    //std::memcpy(&subamp[0], &amplitude[i][int(N/2) - 20], 40 * sizeof(double));
//                    //write(subamp, file);
//                    write(amplitude[i], file);
//                    //pzero.push_back(amplitude[i][int(N / 2)]);
//
//                }
                if (i % samples == 0)
                {
                    k++;
                    if (k% 100 == 0)
                    {
                        std::cout<<k<<endl;
                    }
                    
                    if (k == 2400)
                    {
                        std::cout<<"Hi"<<endl;
                        //std::memcpy(&subamp[0], &amplitude[i][int(N/2) - 20], 40 * sizeof(double));
                        //write(amplitude[i], file);
/*                         for (int l = 0; l < N; l++)
                        {
                            amp.push_back(norm(complexArray[l]));
                        }  */
                        //write(amp, file4);
                        writecomp(file4, complexArray);
                        //cout<<energy[i]<<endl;
                    } 
                    //std::memcpy(&subamp[0], &amplitude[i][int(N/2) - 20], 40 * sizeof(double));
                    //write(subamp, file);
                    //write(amplitude[i], file);

                    pzero.push_back(norm(complexArray[int(N / 2)]));
                    energy.push_back(H(complexArray, alpha, N).real());
                    newnorm.push_back(normm);

                }

            }
        }


        //for (int j = 0; j < steps; j++)
        //{
            //std::memcpy(&subamp[0], &amplitude[j][int(N/2) - 20], 40 * sizeof(double));
            //write(subamp, file);
            //pzero.push_back(amplitude[j][int(N/2)]);
            //if (j % samples == 0)
            //{
            //    std::memcpy(&subamp[0], &amplitude[j][int(N/2) - 20], 40 * sizeof(double));
            //    write(subamp, file);
            //    //pzero.push_back(amplitude[j][int(N / 2)]);
            //}
        //}
        //write(subamp, file);
        writeline(energy, file1);
        write(newnorm, file3);
        //ave.push_back(accumulate( pzero.begin(), pzero.end(), 0.0)/pzero.size());  
        std::cout<<pzero[0]<<endl;            
        write(pzero,file2);
        
        ofstream myfile("./Damping_alpha(" + string(argv[5]) + ")" + "/xi(" + argv[9] + ")" + "/log.txt",std::ofstream::trunc); 
        for (int l = 0; l < 13; l++)
        {
            myfile << argv[l] << endl;
        }
        myfile.close();
    }

    if (model == "Adiabatic")
    {
        string folder = "./Adia";
        struct stat info;
        if (stat(folder.c_str(), &info) != 0)
        {
            mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            string secfolder = folder + "/xi(" + argv[9] + ")";
            mkdir(secfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            file = secfolder + "/lambda(" + argv[5] + ").bin";
            file1 = secfolder + "/energy.bin";
            file2 = secfolder + "/pzero.bin";
            file3 = secfolder + "/norm.bin";
            file4 = secfolder + "/startvec.bin";
	        file5 = secfolder + "/alphavec.txt";
        }
        else if (info.st_mode & S_IFDIR)
        {
            string secfolder = folder + "/xi(" + argv[9] + ")";
            if (stat(secfolder.c_str(), &info) != 0)
            {
                mkdir(secfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
                file = secfolder + "/lambda(" + argv[5] + ").bin";
                file1 = secfolder + "/energy.bin";
                file2 = secfolder + "/pzero.bin";
                file3 = secfolder + "/norm.bin";
                file4 = secfolder + "/startvec.bin";
		        file5 = secfolder + "/alphavec.txt";
                ofstream outFile(file, ios::binary | ios::trunc);
                ofstream outFile1(file1, ios::binary | ios::trunc);
                ofstream outFile2(file2, ios::binary | ios::trunc);
                ofstream outFile3(file3, ios::binary | ios::trunc);
		        ofstream outFile5(file5, ios::binary | ios::trunc);
            }
            else
            {
                file = secfolder + "/lambda(" + argv[5] + ").bin";
                file1 = secfolder + "/energy.bin";
                file2 = secfolder + "/pzero.bin";
                file3 = secfolder + "/norm.bin";
                file4 = secfolder + "/startvec.bin";
		        file5 = secfolder + "/alphavec.txt";
                ofstream outFile(file, ios::binary | ios::trunc);
                ofstream outFile1(file1, ios::binary | ios::trunc);
                ofstream outFile2(file2, ios::binary | ios::trunc);
                ofstream outFile3(file3, ios::binary | ios::trunc);
		        ofstream outFile5(file5, ios::binary | ios::trunc);
            }
        }

        ofstream myfile("./Adia/xi(" + string(argv[9]) + ")" + "/log.txt",std::ofstream::trunc); 

        for (int l = 0; l < 15; l++)
        {
            myfile << argv[l] << endl;
        }
        myfile.close();

        double start = 7.0;
        double end = 3.0;  
        double period = 1;
        double tolerance = 0.1;

        for (double g = start; g > end; g -= period)
        {
            std::cout<<g<<endl;
            for (size_t s = 0; s < ens; s++)
            {
                int k = 0;
                readvec(file4, complexArray, N);
                //cout << s << endl;
                //energy.push_back(H(complexArray, g, N).real());
                energy.push_back(2.0-(g/2));
                newnorm.push_back(c_array_norm(complexArray, N));
                float Norm = newnorm[0];
                for (size_t i = 0; i < steps - 1; i++)
                {  
                    complexArray = EnoughtpartDelta(complexArray, alpha, energy[0], Norm, xi, lambda, dt, N, delta, bond, i);
                    Norm = c_array_norm(complexArray, N);
                    if (i % samples == 0)
                    {  
                        k++;
                        //if (k% 1000 == 0)
                        //{
                        //    cout<<k<<endl;
                        //}

                        if (k == 2400)
                        {
                            writecomp(file4, complexArray);
                        } 

                        if (std::abs(g - 7.99) < tolerance || std::abs(g - 7) < tolerance || std::abs(g - 6) < tolerance || std::abs(g - 5) < tolerance || std::abs(g - 4) < tolerance)
                        {

                            
                            if (k == 2400)
                            {
                                std::cout<<"h1"<<endl;
                                saveComplexVector(file5, complexArray);
                            }
                            pzero.push_back(norm(complexArray[int(N / 2)]));
                            energy.push_back(H(complexArray, g, N).real());
                            newnorm.push_back(Norm);
                        }
                    }
                }
            }
            if (std::abs(g - 7.99) < tolerance || std::abs(g - 7) < tolerance || std::abs(g - 6) < tolerance || std::abs(g - 5) < tolerance || std::abs(g - 4) < tolerance)
            {
                std::cout<<"h2"<<endl;
                writeline(energy, file1);
                writeline(newnorm, file3);
                writeline(pzero,file2);
            }
            energy.clear();
            newnorm.clear();
            pzero.clear();
        }
        //const std::string filepath = "./Adia-0.1/xi(0.8)/alphavec.txt";
        //test = readComplexVector(filepath, complexArray, N);
        
    }
    
    auto end = steady_clock::now();
    auto duration = duration_cast<seconds>(end - start).count();
    std::cout << "duration"
              << "\t" << duration << endl;

    return 0;
}

// defining the Localized wave function
void GaussianArray(vector<complex<double>> &complexArray, double mean, double sigma)
{
    int size = complexArray.size();
    for (int i = 0; i < size; i++)
    {
        double position = i - size / 2;
        double x = position - mean;
        complexArray[i] = exp(-(x * x) / (2 * sigma * sigma));
    }
}

// define a fully localized array
void localizedArray(vector<complex<double>> &complexArray)
{
    int size = complexArray.size();
    for (int i = 0; i < size; i++)
    {
        if (i == int(size / 2))
        {
            complexArray[i] = 1.0;
        }
        else
        {
            complexArray[i] = 0.0;
        }
    }
}

// norm of an array
double c_array_norm(vector<complex<double>> &a, int N)
{
    double result = 0.;
    for (size_t i = 0; i < N; ++i)
    {
        result += norm(a[i]);
    }
    return result;
}

// part of the feedback model
vector<complex<double>> HDNLS(vector<complex<double>> &a, double alpha, double dt, int N)
{
    vector<complex<double>> result(N);

    for (size_t i = 0; i < N; i++)
    {
        int L = i == 0 ? N - 1 : i - 1;
        int R = i == N - 1 ? 0 : i + 1;
        result[i] = I * dt * ((a[R] + a[L] - 2.0 * a[i]) + (alpha * norm(a[i]) * a[i])) + a[i];
    }
//    for (size_t i = 0; i < N; i++)
//    {
//        cout<<result[i]<<" ";
//    }
//    cout<<endl;
    return result;
}

// first damping part
vector<complex<double>> Enoughtpart(vector<complex<double>> &a, double alpha, double Enought, double newnorm, double xi, double lambda, double dt, int N)
{
    vector<complex<double>> result(N);
    complex<double> ham = H(a, alpha, N);
    double common_term = -lambda * (newnorm - 1);

    for (size_t i = 0; i < N; i++)
    {
        int L = (i + N - 1) % N;
        int R = (i + 1) % N;

        result[i] = (I - xi * (ham - Enought)) * dt * ((-a[R] - a[L] + 2.0 * a[i]) - (alpha * norm(a[i]) * a[i])) + a[i] + common_term * a[i];
        

    }

    return result;
}


vector<complex<double>> EnoughtpartDelta(vector<complex<double>> &a, double alpha, double Enought, double newnorm, double xi, double lambda, double dt, int N, double delta, int bond, int place)
{
    vector<complex<double>> result(N);
    complex<double> ham = H(a, alpha, N);
    double common_term = -lambda * (newnorm - 1);

    for (size_t i = 0; i < N; i++)
    {
        int L = (i + N - 1) % N;
        int R = (i + 1) % N;

        if (i < bond || i > (N - bond + 1))
        {
        
            result[i] = (I - xi * (ham - Enought)) * dt * ((-a[R] - a[L] + 2.0 * a[i]) - (alpha * norm(a[i]) * a[i])) - dt * delta * a[i] + a[i] + common_term * a[i];
            
        }
        else
        {
            result[i] = (I - xi * (ham - Enought)) * dt * ((-a[R] - a[L] + 2.0 * a[i]) - (alpha * norm(a[i]) * a[i])) + a[i] + common_term * a[i];
        }

    }

    return result;
}

complex<double> H(vector<complex<double>> &a, double alpha, int N)
{
    complex<double> result = (0,0);

    for (size_t i = 0; i < N; i++)
    {
        int R = i == N - 1 ? 0 : i + 1;
        result = result + (2.0 * norm(a[i]) - a[i] * conj(a[R]) - a[R] * conj(a[i]) - (alpha / 2) * norm(a[i]) * norm(a[i]));
    }

    return result;
}

// adding the noise
void noise(vector<complex<double>> &a, vector<complex<double>> &b, boost::random::mt19937 &rng, boost::normal_distribution<> &m_dist, double gamma, double beta, double dt, int N)
{
    double rand_array1[N];
    double rand_array2[N];
    for (size_t l = 0; l < N; ++l)
    {
        rand_array1[l] = m_dist(rng);
        rand_array2[l] = m_dist(rng);
    }
    for (size_t i = 0; i < N; ++i)
    {
        b[i] = (-I * sqrt(gamma) * (1.0 - beta * norm(a[i])) * (rand_array1[i] + I * rand_array2[i]) * a[i]) - (0.5 * gamma * (1.0 - beta * norm(a[i])) * (1.0 - beta * norm(a[i])) * a[i] * dt);
        // cout<<b[i]<<" ";
    }
    // cout<<endl;
}

// definition of normalization function
void normalize(vector<complex<double>> &complexArray, int N)
{
    double sum = c_array_norm(complexArray, N);
    for (int i = 0; i < N; i++)
    {
        complexArray[i] = complexArray[i] / sqrt(sum);
    }
}

// definition of norm vector function
void normarray(vector<vector<double>> &a, vector<complex<double>> b, int row, double ens)
{
    int size = a[row].size();
    for (int k = 0; k < size; k++)
    {
        a[row][k] = a[row][k] + norm(b[k]) / ens;
    }
}

//void write(vector<double> &A, string &file)
//{
//    std::ofstream outFile(file, std::ios::binary | std::ios::app);
//    outFile.write((char *)A.data(), A.size() * sizeof(double));
//    outFile.close();
//}


void write(const vector<double> &data, const string &file) {
    ofstream outFile(file, ios::binary | ios::trunc);
    for (double value : data) {
        outFile.write(reinterpret_cast<const char*>(&value), sizeof(double));
    }
    outFile.close();
}

void writeline(const vector<double> &data, const string &file) {
    ofstream outFile(file, ios::binary | ios::app); // Use "ios::app" to append to the file
    for (double value : data) {
        outFile.write(reinterpret_cast<const char*>(&value), sizeof(double));
    }
    outFile.close();
}

void writecomp(const std::string& file_path, vector<complex<double>>& complex_vector) {
    std::ofstream fille(file_path, std::ios::binary| std::ios::trunc);

    if (!fille) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return;
    }

    // Write each complex value to the binary file
    for (const auto& complex_value : complex_vector) {
        fille.write(reinterpret_cast<const char*>(&complex_value), sizeof(std::complex<double>));
    }
}

void writecompadd(const std::string& file_path, vector<complex<double>>& complex_vector) {
    std::ofstream fille(file_path, std::ios::binary| std::ios::app);

    if (!fille) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return;
    }

    // Write each complex value to the binary file
    for (const auto& complex_value : complex_vector) {
        fille.write(reinterpret_cast<const char*>(&complex_value), sizeof(std::complex<double>));
    }
    fille.close();
}

void readvec(const std::string& filePath, vector<complex<double>>& complexVector, size_t N) {
    // Open the binary file for reading in binary mode
    ifstream file(filePath, ios::binary);

    // Check if the file is opened successfully
    if (!file.is_open()) {
        cerr << "Error: Unable to open file." << endl;
        return;
    }

    // Read the complex values from the file and fill the vector
    for (size_t i = 0; i < N; ++i) {
        double realPart, imagPart;

        // Read the real and imaginary parts from the file
        if (!file.read(reinterpret_cast<char*>(&realPart), sizeof(double))) {
            cerr << "Error reading real part from file." << endl;
            break;
        }

        if (!file.read(reinterpret_cast<char*>(&imagPart), sizeof(double))) {
            cerr << "Error reading imaginary part from file." << endl;
            break;
        }

        // Create a complex number and add it to the vector
        complex<double> complexValue(realPart, imagPart);
        //complexVector.push_back(complexValue);
        complexVector[i] = complexValue;
    }

    // Close the file
    file.close();
}


void NormalizedGaussian(std::vector<std::complex<double>>& result, double peak, double standardDeviation) {
    int size = result.size();
    double sum = 0.0;

    if (size % 2 == 0) {
        std::cerr << "Size should be odd for a centered peak." << std::endl;
        return;
    }

    int middle = size / 2;

    for (int i = 0; i < size; ++i) {
        int x = i - middle;
        double realPart = peak * exp(-0.5 * (x * x) / (standardDeviation * standardDeviation));
        result[i] = std::complex<double>(realPart, 0.0);
        sum += norm(result[i]);
    }

    // Normalize the real part of the array to ensure that the sum of values is equal to 1.
    for (int i = 0; i < size; ++i) {
        result[i] /= sqrt(sum);
    }
}


void saveComplexVector(const std::string& file_path, const std::vector<std::complex<double>>& complex_vector) {
    // Open the file in append mode to add data if it already exists
    std::ofstream file(file_path, std::ios::app);

    if (!file) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return;
    }

    // Write each complex value to the text file
    for (const auto& complex_value : complex_vector) {
        file << complex_value.real() << " " << complex_value.imag() << "\n";
    }

    // Close the file
    file.close();
}

std::vector<std::complex<double>> readComplexVector(const std::string& file_path) {
    std::vector<std::complex<double>> complex_vector;

    // Open the file for reading
    std::ifstream file(file_path);

    if (!file) {
        std::cerr << "Error opening file: " << file_path << std::endl;
        return complex_vector;  // Return an empty vector on error
    }

    double real_part, imag_part;

    // Read data from the file and store it in the complex vector
    while (file >> real_part >> imag_part) {
        //cout << real_part << imag_part << endl;
        complex_vector.emplace_back(real_part, imag_part);
    }

    // Close the file
    file.close();

    return complex_vector;
}


