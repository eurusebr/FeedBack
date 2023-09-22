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
void write(vector<double> &A, string &file);
vector<complex<double>> Enoughtpart(vector<complex<double>> &a, double alpha, double Enought, double newnorm, double xi, double lambda, double dt, int N);
complex<double> H(vector<complex<double>> &a, double alpha, int N);
vector<complex<double>> EnoughtpartDelta(vector<complex<double>> &a, double alpha, double Enought, double newnorm, double xi, double lambda, double dt, int N, double delta, int bond);

int main(int argc, char **argv)
{

#include <iostream>
#include <vector>
#include <sstream>

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
    
    string file, file1, file2;
    vector<vector<double>> amplitude(steps, vector<double>(N, 0));
    vector<double> pzero, energy, ave;
    vector<double> subamp(41);
    vector<double> subbamp(N);

    boost::random::mt19937 rng;
    boost::normal_distribution<> m_dist(0.0, sqrt(0.5 * dt));
    vector<complex<double>> complexArray(N);
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
        // Create a new directory
        string folder = "./DNLS";
        struct stat info;

        if (stat(folder.c_str(), &info) != 0)
        {
            mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            file = folder + "/" + argv[5] + ".bin";
            file1 = folder + "/energy.bin";
            file2 = folder + "/pzero(" + argv[5] + ").bin";
        }

        if (stat(folder.c_str(), &info) == 0)
        {
            file = folder + "/" + argv[5] + ".bin";
            ofstream outFile(file, ios::binary | ios::trunc);
            file1 = folder + "/energy.bin";
            file2 = folder + "/pzero(" + argv[5] + ").bin";
        }

        // GaussianArray(complexArray, 0.0, sigma);
        localizedArray(complexArray);
        normalize(complexArray, N);
        //energy.push_back(H(complexArray, alpha, N).real());
        pzero.push_back(amplitude[0][int(N / 2)]);
        normarray(amplitude, complexArray, 0, ens);
        //std::memcpy(&subamp[0], &amplitude[0][int(N / 2) - 20], 40 * sizeof(double));
        //write(subamp, file);

        for (size_t j = 0; j < steps-1; j++)
        {
            complexArray = HDNLS(complexArray, alpha, dt, N);
            normalize(complexArray, N);
            
            normarray(amplitude, complexArray, j+1, ens);
//            write(subamp, file);
            if (j % samples == 0)
            {
                energy.push_back(H(complexArray, alpha, N).real());
                pzero.push_back(amplitude[j][int(N / 2)]);
                std::memcpy(&subamp[0], &amplitude[j][int(N / 2) - 20], 40 * sizeof(double));
                write(subamp, file);
            }
        }
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
            cout << s << endl;
            localizedArray(complexArray);
            normalize(complexArray, N);
            normarray(amplitude, complexArray, 0, ens);
            for (size_t i = 0; i < steps - 1; i++)
            {
                noise(complexArray, AddingArray, rng, m_dist, gamma, beta, dt, N);
                complexArray = HDNLS(complexArray, alpha, dt, N);
                for (size_t j = 0; j < N; j++)
                {
                    complexArray[j] = complexArray[j] + AddingArray[j];
                }
                normalize(complexArray, N);
                normarray(amplitude, complexArray, i + 1, ens);
            }
        }

        cout << steps << endl;
        for (int j = 0; j < steps; j++)
        {
            // std::memcpy(&subamp[0], &amplitude[j][int(N/2) - 20], 40 * sizeof(double));
            // pzero.push_back(amplitude[j][int(N/2)]);
            // write(subamp, file);
            if (j % samples == 0)
            {
                pzero.push_back(amplitude[j][int(N / 2)]);
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
                ofstream outFile(file, ios::binary | ios::trunc);
                ofstream outFile1(file1, ios::binary | ios::trunc);
                ofstream outFile2(file2, ios::binary | ios::trunc);
            }
            else
            {
                file = secfolder + "/lambda(" + argv[8] + ").bin";
                file1 = secfolder + "/energy.bin";
                file2 = secfolder + "/pzero.bin";
                ofstream outFile(file, ios::binary | ios::trunc);
                ofstream outFile1(file1, ios::binary | ios::trunc);
                ofstream outFile2(file2, ios::binary | ios::trunc);
            }
        }


        for (size_t s = 0; s < ens; s++)
        {
            cout << s << endl;
            localizedArray(complexArray);
            double newnorm = c_array_norm(complexArray, N);
            //normalize(complexArray, N);
            energy.push_back(H(complexArray, alpha, N).real());
            normarray(amplitude, complexArray, 0, ens);
            //std::memcpy(&subamp[0], &amplitude[0][int(N/2) - 20], 40 * sizeof(double));
            //write(subamp, file);
            for (size_t i = 0; i < steps - 1; i++)
            {
                ///complexArray = Enoughtpart(complexArray, alpha, energy[0], newnorm, xi, lambda, dt, N);
                complexArray = EnoughtpartDelta(complexArray, alpha, energy[0], newnorm, xi, lambda, dt, N, delta, bond);
                newnorm = c_array_norm(complexArray, N);
                //normalize(complexArray, N);
                //energy.push_back(H(complexArray, alpha, N).real());
                normarray(amplitude, complexArray, i + 1, ens);
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
                    //std::memcpy(&subamp[0], &amplitude[i][int(N/2) - 20], 40 * sizeof(double));
                    //write(subamp, file);
                    //write(amplitude[i], file);
                    pzero.push_back(amplitude[i][int(N / 2)]);

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
        //write(energy, file1);
        //ave.push_back(accumulate( pzero.begin(), pzero.end(), 0.0)/pzero.size());              
        write(pzero,file2);
        
        ofstream myfile("./Damping_alpha(" + string(argv[5]) + ")" + "/xi(" + argv[9] + ")" + "/log.txt",std::ofstream::trunc); 
        for (int l = 0; l < 13; l++)
        {
            myfile << argv[l] << endl;
        }
        myfile.close();
    }

    if (model == "DampingAlpha")
    {
       vector<float> myVector{ 1, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.50, 3.55, 3.60, 3.65, 3.70, 3.75, 3.80, 3.85, 3.90, 3.95, 4, 4.25, 4.5, 4.75, 5, 5.25, 5.5, 5.75, 6, 6.25, 6.5, 6.75, 7, 7.25, 7.5, 7.75, 8 }; 
       
        // Create a new directory
        string folder = "./AllAlpha";
        struct stat info;
        if (stat(folder.c_str(), &info) != 0)
        {
            mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            file = folder + "/ave.bin";
            file1 = folder + "/dens.bin";
            file2 = folder + "/energy.bin";


        }
        else if (info.st_mode & S_IFDIR)
        {
            file = folder + "/ave.bin";
            ofstream outFile(file, ios::binary | ios::trunc);
            file1 = folder + "/dens.bin";
            ofstream outFile1(file1, ios::binary | ios::trunc);
            file2 = folder + "/energy.bin";
            ofstream outFile2(file2, ios::binary | ios::trunc);
        }

        for (size_t s = 0; s < myVector.size(); s++)
        {
            cout << s << endl;
            localizedArray(complexArray);
            double newnorm = c_array_norm(complexArray, N);
            energy.push_back(H(complexArray, myVector[s], N).real());
            normarray(amplitude, complexArray, 0, ens);
            for (size_t i = 0; i < steps - 1; i++)
            {
                complexArray = Enoughtpart(complexArray, myVector[s], energy[0], newnorm, xi, lambda, dt, N);
                newnorm = c_array_norm(complexArray, N);
                normarray(amplitude, complexArray, i + 1, ens);
                if (i % samples == 0)
                {
                    //std::memcpy(&subamp[0], &amplitude[i][int(N/2) - 20], 40 * sizeof(double));
                    //write(subamp, file);
                    pzero.push_back(amplitude[i][int(N / 2)]);
                    //write(amplitude[i], file);
                }
                //if (i == 4000000)
                //{
                //    std::memcpy(&subbamp[0], &amplitude[i][int(N/2) - 250], N * sizeof(double));
                //    write(subbamp, file1);
                    //pzero.push_back(amplitude[i][int(N / 2)]);
                //}

            }
            write(energy,file2);
            //ofstream pzz("./AllAlpha/pzero" + std::to_string(s) + ".txt",std::ofstream::trunc); 
            //for (int l = 0; l < pzero.size(); l++)
            //{
            //    pzz << pzero[l] << endl;
            //}
            //pzz.close();
            ave.push_back(accumulate( pzero.begin(), pzero.end(), 0.0)/pzero.size());   
            pzero.clear();
            energy.clear();
            for(auto& elem : amplitude) std::fill(elem.begin(), elem.end(), 0);        
        }

//        for (size_t i = 0; i < 2; i++)
//        {
//            cout<<ave[i]<<' ';
//        }
//        cout<<endl;      

        write(ave,file);

        ofstream myfile("./AllAlpha/log.txt",std::ofstream::trunc); 
        for (int l = 0; l < 13; l++)
        {
            myfile << argv[l] << endl;
        }
        myfile.close();
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

    for (size_t i = 0; i < N; i++)
    {
        int L = i == 0 ? N - 1 : i - 1;
        int R = i == N - 1 ? 0 : i + 1;
        result[i] = (I - xi * (ham - Enought)) * dt * ((-a[R] - a[L] + 2.0 * a[i]) - (alpha * norm(a[i]) * a[i])) - lambda * a[i] * (newnorm - 1) + a[i];
        
    }
    return result;
}

vector<complex<double>> EnoughtpartDelta(vector<complex<double>> &a, double alpha, double Enought, double newnorm, double xi, double lambda, double dt, int N, double delta, int bond)
{
    vector<complex<double>> result(N);
    complex<double> ham = H(a, alpha, N);    

    for (size_t i = 0; i < N; i++)
    {
        if(i < bond || i > (N-bond)) 
        {
            int L = i == 0 ? N - 1 : i - 1;
            int R = i == N - 1 ? 0 : i + 1;
            result[i] = (I - xi * (ham - Enought)) * dt * ((-a[R] - a[L] + 2.0 * a[i]) - (alpha * norm(a[i]) * a[i])) - lambda * a[i] * (newnorm - 1) - (I * delta * a[i]) + a[i];
        }
        else
        {
            int L = i == 0 ? N - 1 : i - 1;
            int R = i == N - 1 ? 0 : i + 1;
            result[i] = (I - xi * (ham - Enought)) * dt * ((-a[R] - a[L] + 2.0 * a[i]) - (alpha * norm(a[i]) * a[i])) - lambda * a[i] * (newnorm - 1) + a[i];
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

void write(vector<double> &A, string &file)
{
    std::ofstream outFile(file, std::ios::binary | std::ios::app);
    outFile.write((char *)A.data(), A.size() * sizeof(double));
    outFile.close();
}
