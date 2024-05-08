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
#include <sstream>
using namespace std;
using namespace boost::numeric::odeint;
using namespace std::chrono;

// definition of imaginary unit
const complex<double> I(0.0, 1.0);

vector<complex<double>> EnoughtpartDelta(vector<complex<double>> &a, double alpha, double Enought, double newnorm, double xi, double lambda, double dt, int N, double delta, int bond, int place);
void writeline(const vector<double> &data, const string &file);
void writecomp(const std::string& file_path, vector<complex<double>>& complex_vector);
void write(const vector<double> &data, const string &file);
complex<double> H(vector<complex<double>> &a, double alpha, int N);
double c_array_norm(vector<complex<double>> &a, int N);
void localizedArray(vector<complex<double>> &complexArray);
vector<complex<double>> Enoughtpart(vector<complex<double>> &a, double alpha, double Enought, double newnorm, double xi, double lambda, double dt, int N);
void readvec(const std::string& filePath, vector<complex<double>>& complexVector, size_t N);
void saveComplexVector(const std::string& file_path, const std::vector<std::complex<double>>& complex_vector);

int main(int argc, char **argv)
{

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
    double xi = stod(argv[9]);
    int steps = stod(argv[10]);
    double ens = stod(argv[11]);
    int samples = stod(argv[12]);
    double delta = stod(argv[13]); 
    int bond = stod(argv[14]);

    auto start = steady_clock::now();
    
    string file, file1, file2, file3, file4, file5, file6, file7;
    vector<double> pzero, energy, pone, pminusone;
    vector<double> newnorm;
    vector<complex<double>> complexArray(N);
    vector<complex<double>> midcompvec;

// Create a new directory
    string folder = "./Omega";
    struct stat info;
    if (stat(folder.c_str(), &info) != 0)
    {
        mkdir(folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        string secfolder = folder + "/alpha(" + argv[5] + ")";
        mkdir(secfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        file = secfolder + "/alphavec.txt";
        file1 = secfolder + "/energy.bin";
        file2 = secfolder + "/pzero.bin";
        file3 = secfolder + "/norm.bin";
        file4 = secfolder + "/startvec.txt";
        file5 = secfolder + "/pone.bin";
        file6 = secfolder + "/pminusone.bin";
        file7 = secfolder + "/compmid.bin";
    }
    else if (info.st_mode & S_IFDIR)
    {
        string secfolder = folder + "/alpha(" + argv[5] + ")";
        if (stat(secfolder.c_str(), &info) != 0)
        {
            mkdir(secfolder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
            file = secfolder + "/alphavec.txt";
            file1 = secfolder + "/energy.bin";
            file2 = secfolder + "/pzero.bin";
            file3 = secfolder + "/norm.bin";
            file4 = secfolder + "/startvec.txt";
            file5 = secfolder + "/pone.bin";
            file6 = secfolder + "/pminusone.bin";
            file7 = secfolder + "/compmid.bin";
            ofstream outFile(file, ios::binary | ios::trunc);
            ofstream outFile1(file1, ios::binary | ios::trunc);
            ofstream outFile2(file2, ios::binary | ios::trunc);
            ofstream outFile3(file3, ios::binary | ios::trunc);
            ofstream outFile5(file5, ios::binary | ios::trunc);
            ofstream outFile6(file6, ios::binary | ios::trunc);
            ofstream outFile7(file7, ios::binary | ios::trunc);
        }
        else
        {
            file = secfolder + "/alphavec.txt";
            file1 = secfolder + "/energy.bin";
            file2 = secfolder + "/pzero.bin";
            file3 = secfolder + "/norm.bin";
            file4 = secfolder + "/startvec.txt";
            file5 = secfolder + "/pone.bin";
            file6 = secfolder + "/pminusone.bin";
            file7 = secfolder + "/compmid.bin";
            ofstream outFile(file, ios::binary | ios::trunc);
            ofstream outFile1(file1, ios::binary | ios::trunc);
            ofstream outFile2(file2, ios::binary | ios::trunc);
            ofstream outFile3(file3, ios::binary | ios::trunc);
            ofstream outFile5(file5, ios::binary | ios::trunc);
            ofstream outFile6(file6, ios::binary | ios::trunc);
            ofstream outFile7(file7, ios::binary | ios::trunc);
        }
    }

    std::string filename = std::string("./Omega/alpha(") + argv[5] + ")/log.txt";
    std::ofstream myfile(filename.c_str(), std::ofstream::trunc);
    for (int l = 0; l < 15; l++)
    {
        myfile << argv[l] << endl;
    }
    myfile.close();
    std::vector<double> Enot = {-1.78878, -1.33827, -0.908747, -0.518153};
    int count = 0;
    readvec(file4, complexArray, N);
    newnorm.push_back(c_array_norm(complexArray, N));
    double normm = newnorm[0];
    energy.push_back(Enot[count]);
    for (size_t s = 0; s < ens; s++)
    {
        //std::cout << s << endl;
        int k = 0;
        for (size_t i = 0; i < steps - 1; i++)
        {
            complexArray = EnoughtpartDelta(complexArray, alpha, energy[0], normm, xi, lambda, dt, N, delta, bond, i);
            normm = c_array_norm(complexArray, N);
            if (i % samples == 0)
            {
                k++;
                //if (k% 100 == 0)
                //{
                //    std::cout<<k<<endl;
                //}
                
                //if (k == 2400)
                //{
                //    cout<<"Hi"<<endl;
                //    writecomp(file4, complexArray);
                //} 
                pzero.push_back(norm(complexArray[int(N / 2)]));
                pone.push_back(norm(complexArray[int(N / 2) + 1]));
                pminusone.push_back(norm(complexArray[int(N / 2) - 1]));
                energy.push_back(H(complexArray, alpha, N).real());
                newnorm.push_back(normm);
                midcompvec.push_back(complexArray[int(N / 2)]);
            }
        }
    }
    cout<< norm(complexArray[int(N / 2)]) << " " << complexArray[int(N / 2)] <<endl;
    saveComplexVector(file, complexArray);
    saveComplexVector(file7, midcompvec);
    cout<<"h2"<<endl;
    writeline(energy, file1);
    writeline(newnorm, file3);
    std::cout<<"pzero: "<<pzero[0]<<endl;    
    std::cout<<"energy: "<<energy.back()<<endl;                    
    writeline(pzero,file2);
    writeline(pone,file5);
    writeline(pminusone,file6);  
    writecomp(file4, complexArray);               
    
    auto end = steady_clock::now();
    auto duration = duration_cast<seconds>(end - start).count();
    std::cout << "duration"
              << "\t" << duration << endl;

    return 0;

}

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

void write(const vector<double> &data, const string &file) {
    ofstream outFile(file, ios::binary | ios::trunc);
    for (double value : data) {
        outFile.write(reinterpret_cast<const char*>(&value), sizeof(double));
    }
    outFile.close();
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

double c_array_norm(vector<complex<double>> &a, int N)
{
    double result = 0.;
    for (size_t i = 0; i < N; ++i)
    {
        result += norm(a[i]);
    }
    return result;
}

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

void readvec(const std::string& filePath, vector<complex<double>>& complexVector, size_t N) {
    // Open the text file for reading
    ifstream file(filePath);

    // Check if the file is opened successfully
    if (!file.is_open()) {
        cerr << "Error: Unable to open file." << endl;
        return;
    }

    // Read the complex values from the file and fill the vector
    for (size_t i = 0; i < N; ++i) {
        double realPart, imagPart;

        // Read the real and imaginary parts from the file
        if (!(file >> realPart >> imagPart)) {
            cerr << "Error reading complex number from file." << endl;
            break;
        }

        // Create a complex number and add it to the vector
        complex<double> complexValue(realPart, imagPart);
        complexVector[i] = complexValue;
    }
    
    // Close the file
    file.close();
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