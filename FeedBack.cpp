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

//definition of imaginary unit
const complex <double> I(0.0, 1.0);


// defining functions
void GaussianArray(vector<complex<double>>& complexArray, double mean, double sigma);
double c_array_norm(vector<complex<double>>& a, int N);
vector<complex<double>> HDNLS(vector<complex<double>>& a, double alpha, double dt, int N);
void noise(vector<complex<double>>& a, vector<complex<double>>& b, boost::random::mt19937 & rng, boost::normal_distribution<> & m_dist, double gamma, double beta, double dt, int N);
void normalize(vector<complex<double>> &complexArray, int N);
void localizedArray(vector<complex<double>>& complexArray);
void normarray(vector<vector<double>>& a, vector<complex<double>> b, int row, double ens);



int main(int argc, char** argv)
{
    auto start = steady_clock::now();   
    string model = argv[1];
    int N = stod(argv[2]);
    double sigma = stod(argv[3]);
    double dt = stod(argv[4]);
    double alpha = stod(argv[5]);
    double gamma = stod(argv[6]);
    double beta = stod(argv[7]);
    int steps = stod(argv[8]);
    double ens = stod(argv[9]);
    string file;
    vector<vector<double>> amplitude(steps,vector<double>(N,0));
    vector<double> pzero;
    vector<double> subamp(41);

    boost::random::mt19937 rng;                                        
    boost::normal_distribution<> m_dist(0.0, sqrt(0.5*dt));
    vector<complex<double>> complexArray(N);
    vector<complex<double>> noiseArray(N);
    vector<complex<double>> AddingArray(N);

    // GaussianArray(complexArray, 0.0, sigma);
    // cout<<"complex:"<<endl;
    // for (size_t s = 0; s < N; s++)
    // {
    //     std::cout << complexArray[s]<<" ";
    // }
    // std::cout << endl;
    // pzero.push_back(amplitude[int(N/2)]);

    if(model == "DNLS")
    {
        // Create a new directory 
        string folder = "./DNLS";
        struct stat info;
    
        if( stat( folder.c_str(), &info ) != 0 ) 
        {
            mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            file = folder + "/" + argv[5] + ".bin";
        }

        if( stat( folder.c_str(), &info ) == 0 )
        {
            file = folder + "/" + argv[5] + ".bin";
            ofstream outFile(file, ios::binary | ios::trunc);
        }

        localizedArray(complexArray);
        normalize(complexArray, N);
        normarray(amplitude, complexArray, 0, ens);
        std::memcpy(&subamp[0], &amplitude[0][int(N/2) - 20], 40 * sizeof(double));
        std::ofstream outFile(file, std::ios::binary | std::ios::app);
        outFile.write((char*)subamp.data(), subamp.size() * sizeof(double));
        outFile.close();

        for (size_t j = 0; j < steps; j++)
        {
            complexArray = HDNLS(complexArray, alpha, dt, N);
            normalize(complexArray, N);
            normarray(amplitude, complexArray, j, ens);
            std::memcpy(&subamp[0], &amplitude[int(N/2) - 20], 40 * sizeof(double));

            std::ofstream outFile(file, std::ios::binary | std::ios::app);
            outFile.write((char*)subamp.data(), subamp.size() * sizeof(double));
            outFile.close();
        }

        // std::ofstream outFile(file, std::ios::binary | std::ios::app);
        // outFile.write((char*)pzero.data(), pzero.size() * sizeof(double));
        // outFile.close();
    }

    if(model == "FeedBack")
    {
        // Create a new directory 
        string folder = "./FeedBack_beta(" + string(argv[7]) + ")";
        struct stat info;
        if( stat( folder.c_str(), &info ) != 0 )
        {
            mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            string secfolder = folder + "/alpha(" + argv[5] + ")";
            mkdir(secfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            file = secfolder + "/gamma(" + argv[6] + ").bin";
        }
        else if( info.st_mode & S_IFDIR )
        {
            string secfolder = folder + "/alpha(" + argv[5] + ")";
            if( stat( secfolder.c_str(), &info ) != 0 )
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
            cout<<s<<endl;
            localizedArray(complexArray);
            normalize(complexArray, N);
            normarray(amplitude, complexArray, 0, ens);
            for (size_t i = 0; i < steps-1; i++)
            {
                noise(complexArray, AddingArray, rng, m_dist, gamma, beta, dt, N);
                complexArray = HDNLS(complexArray, alpha, dt, N);
                for (size_t j = 0; j < N; j++)
                {
                    complexArray[j] = complexArray[j] + AddingArray[j];
                }
                normalize(complexArray, N);
                normarray(amplitude, complexArray, i+1, ens);
            }
        }
        
        // std::ofstream outFile(file, std::ios::binary | std::ios::app);
        for (int i = 0; i < steps; i++) 
        {
            std::memcpy(&subamp[0], &amplitude[i][int(N/2) - 20], 40 * sizeof(double));
            pzero.push_back(amplitude[i][int(N/2)]);
            // outFile.write((char*)subamp.data(), subamp.size() * sizeof(double));
        }
        // outFile.close();

        std::ofstream outFile(file, std::ios::binary | std::ios::app);
        outFile.write((char*)pzero.data(), pzero.size() * sizeof(double));
        outFile.close();

    }

    // std::ofstream outFile(file, std::ios::binary | std::ios::app);
    // outFile.write((char*)pzero.data(), pzero.size() * sizeof(double));
    // outFile.close();

    auto end = steady_clock::now();
    auto duration = duration_cast<seconds>(end - start).count();
    std::cout << "duration" << "\t" << duration << endl;
    
    return 0;

}

// defining the Localized wave function
void GaussianArray(vector<complex<double>>& complexArray, double mean, double sigma)
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
void localizedArray(vector<complex<double>>& complexArray)
{
    int size = complexArray.size();
    for (int i = 0; i < size; i++) 
    {
        if(i == int(size/2))
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
double c_array_norm(vector<complex<double>>& a, int N) 
{
    double result = 0.;
    for (size_t i = 0; i < N; ++i) 
    {
        result += norm(a[i]);
    }
    return result;
}

// part of the feedback model
vector<complex<double>> HDNLS(vector<complex<double>>& a, double alpha, double dt, int N) 
{
    vector<complex<double>> result(N);

    for (size_t i = 0; i < N; i++)
    {
        int L = i==0 ? N-1 : i-1;
        int R = i==N-1 ? 0 : i+1;
        result[i] = I * dt *((a[R] + a[L] - 2.0 * a[i]) + (alpha * norm(a[i]) * a[i])) + a[i];
    }
    // cout<<endl;
    return result;
}


// adding the noise
void noise(vector<complex<double>>& a, vector<complex<double>>& b, boost::random::mt19937 & rng, boost::normal_distribution<> & m_dist, double gamma, double beta, double dt, int N)
{
    double rand_array1[N];
    double rand_array2[N];
    for(size_t l = 0; l < N; ++l) 
    {
        rand_array1[l] = m_dist(rng);
        rand_array2[l] = m_dist(rng);
    }
    for (size_t i = 0; i < N; ++i) 
    {
        b[i] = (-I * sqrt(gamma) * (1.0 - beta * norm(a[i])) * (rand_array1[i] + I * rand_array2[i]) * a[i]) - (0.5 * gamma * (1.0 - beta * norm(a[i]))* (1.0 - beta * norm(a[i])) * a[i] * dt) ;
        // cout<<b[i]<<" ";
    }
    // cout<<endl;
}


//definition of normalization function
void normalize(vector<complex<double>>& complexArray, int N)
{
    double sum = c_array_norm(complexArray, N) ;
    for (int i = 0; i < N; i++)
    {
        complexArray[i] = complexArray[i]/sqrt(sum);
    }
}

//definition of norm vector function
void normarray(vector<vector<double>>& a, vector<complex<double>> b, int row, double ens)
{
    int size = a[row].size();
    for (int k = 0; k < size; k++) 
    {
        a[row][k] = a[row][k] + norm(b[k])/ens;
    }
}
