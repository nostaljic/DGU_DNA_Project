////////////////////////////////////////////////////////
//////                                  /////
//////    동국대학교 컴퓨터공학과             /////
//////    수업 : 컴퓨터알고리즘과실습          /////
//////    학번 : 2016112096                    /////
//////    이름 : 박재용                      /////
//////    <scorpion@dgu.ac.kr>                /////
//////                                  /////
//////    문제 3                           /////
//////                                  /////
////////////////////////////////////////////////////////
//#include <iostream>
//#include <fstream>
//#include <string>
//#include <vector>
//#include <chrono>
//#include <time.h>
//#include <random>
//#include <algorithm>
//#include <random>
//#include <ctime>
//#include <functional>
//
//using namespace std;
//vector<char> TWIG = { 'A','T','G','C' };
//vector<int> randomSeed = { 7,13,17,23 };
////chrono::system_clock::time_point start = chrono::system_clock::now();
////chrono::duration<double> sec = chrono::system_clock::now() - start;
//vector<string> shortreads_container;
///*
//일반적인 랜덤함수를 쓰게 될 경우 문자열(DNA)의 반복 정도가 심해집니다.
//반복을 줄일 수 있는방법으로 TWIG를 int seed 값으로 모듈러 연산을 한 값이 0이 되는 순간 섞어줍니다.
//그 후 seed를 randomSeed 인트 벡터의 숫자를 더해 업데이트 해줍니다.
//*/
////auto make_string_file = [](int fileLength, string filePath) {
////    ofstream dna(filePath);
////    string fileDnaString;
////    random_device rd;
////    int seed = 1007;
////    for (int i = 0; i < fileLength; i++) {
////        fileDnaString += TWIG[(rand() % 4)];
////        if (i % seed == 0) {
////            shuffle(TWIG.begin(), TWIG.end(), default_random_engine(rd()));
////            seed += randomSeed[rand() % 4];
////        }
////    }
////    dna << fileDnaString;
////    dna.close();
////};
//////Reference DNA로부터 서브시퀀스의 길이보다 작은 곳에서의 Twig를 변경해준다.
////auto make_my_dna = [](int subSequenceLength, int short_reads_count, int snp_count, int fileLength, string refDnaPath, string myDnaPath) {
////
////    ifstream refDNA(refDnaPath);
////    string refDnaString;
////    getline(refDNA, refDnaString);
////
////    ofstream myDNA(myDnaPath);
////
////    for (int i = 0; i < fileLength - subSequenceLength - 1; i += subSequenceLength) {
////        for (int j = 0; j < snp_count - 1; j++) {
////            refDnaString[i + (rand() % subSequenceLength)] = TWIG[(rand() % 4)];
////        }
////    }
////    myDNA << refDnaString;
////
////    refDNA.close();
////    myDNA.close();
////};
////////제한 효소로 원본 데이터(Reference DNA sequence)로부터 절단된 short reads를 얻는다.
//////auto endo_nuclease = [](int subSequenceLength, int shorReadsCounter, int fileLength, string filePath, string shortReadsPath) {
//////   //vector<string> shortReads;
//////   ifstream dna(filePath);
//////   string fileDnaString;
//////   getline(dna, fileDnaString);
//////
//////   ofstream shortReads(shortReadsPath);
//////   int num_of_short_reads = 0;
//////   int sign_flag = 1;
//////   int index_of_sub_dna = 0;
//////   //int rand_degree = 1;
//////   int delim = fileLength / shorReadsCounter;
//////   while (num_of_short_reads < shorReadsCounter - 1) {
//////      if (num_of_short_reads != 0)index_of_sub_dna += sign_flag * ((rand() % (subSequenceLength - delim)) + delim);
//////      if (index_of_sub_dna + subSequenceLength > fileLength - 1) { sign_flag = -1;  continue; }
//////      if (index_of_sub_dna < 0) { sign_flag = 1;  continue; }
//////      shortReads << fileDnaString.substr(index_of_sub_dna, subSequenceLength); shortReads << '\n';
//////      shortReads << fileDnaString.substr(fileLength - (subSequenceLength + index_of_sub_dna), subSequenceLength); shortReads << '\n';
//////      shortreads_container.push_back(fileDnaString.substr(index_of_sub_dna, subSequenceLength));
//////      shortreads_container.push_back(fileDnaString.substr(fileLength - (subSequenceLength + index_of_sub_dna), subSequenceLength));
//////      //num_of_short_reads += 1;
//////      num_of_short_reads += 2;
//////   }
//////   shortReads << fileDnaString.substr(fileLength - subSequenceLength, subSequenceLength); shortReads << '\n';
//////   shortreads_container.push_back(fileDnaString.substr(fileLength - subSequenceLength, subSequenceLength));
//////
//////   shortReads.close();
//////   dna.close();
//////};
//////제한 효소로 원본 데이터(Reference DNA sequence)로부터 절단된 short reads를 얻는다.
////auto endo_nuclease = [](int subSequenceLength, int shortReadsCounter, int fileLength, string filePath, string shortReadsPath) {
////    //vector<string> shortReads;
////    ifstream dna(filePath);
////    string fileDnaString;
////    getline(dna, fileDnaString);
////
////    ofstream shortReads(shortReadsPath);
////
////    mt19937 engine((unsigned int)time(NULL));                    // 松本眞のアルゴリズム
////    uniform_int_distribution<int> distribution(0, fileLength - subSequenceLength - 1);       // 생성 범위
////    auto generator = bind(distribution, engine);
////    //int maxind = fileLength - subSequenceLength;
////    for (int i = 0; i < shortReadsCounter; i++) {
////        int randomindex = generator();
////        //int randomindex = rand()% maxind;
////        shortReads << fileDnaString.substr(randomindex, subSequenceLength); shortReads << '\n';
////        shortreads_container.push_back(fileDnaString.substr(randomindex, subSequenceLength));
////    }
////
////    shortReads.close();
////    dna.close();
////};
////auto endo_nuclease = [](int subSequenceLength, int shortReadsCounter, int fileLength, string filePath, string shortReadsPath) {
////   //vector<string> shortReads;
////   ifstream dna(filePath);
////   string fileDnaString;
////   getline(dna, fileDnaString);
////
////   ofstream shortReads(shortReadsPath);
////   srand(time(NULL));
////   //mt19937 engine((unsigned int)time(NULL));                    // 松本眞のアルゴリズム
////   //uniform_int_distribution<int> distribution(0, fileLength - subSequenceLength);       // 생성 범위
////   //auto generator = bind(distribution, engine);
////   int maxind = fileLength - subSequenceLength;
////   for (int i = 0; i < shortReadsCounter; i++) {
////      
////      int randomindex = (rand()*97)%maxind;
////      //cout << randomindex << endl;
////      shortReads << fileDnaString.substr(randomindex, subSequenceLength); shortReads << '\n';
////      shortreads_container.push_back(fileDnaString.substr(randomindex, subSequenceLength));
////   }
////
////   shortReads.close();
////   dna.close();
////};
////brute force 방식으로 복원한다.
//auto trivial_reconstruction = [](int subSequenceLength, int snp_count, int fileLength, string reffilePath, string shortReadsPath, string reconstructedfilePath) {
//    //vector<string> shortReads;
//    chrono::system_clock::time_point start = chrono::system_clock::now();
//    ofstream reconstructed_dna(reconstructedfilePath);
//    string reconstructed; reconstructed.resize(fileLength + 1);
//    //cout << reconstructed << endl;
//    ifstream RefDna(reffilePath);
//    string refDNA;
//    getline(RefDna, refDNA);
//
//    ifstream shortReads(shortReadsPath);
//    string shortRead;
//
//    int count_shortread = 0;
//    //while (getline(shortReads, shortRead)) {
//
//    int delim = fileLength / (shortreads_container.size());
//    while (getline(shortReads, shortRead)) {
//        //shortRead = shortreads_container[count_shortread];
//        if (count_shortread++ % 1000 == 0)cout << (double)count_shortread / 400000.0 << "\n";
//        //cout << shortRead << endl;
//
//        int snp_flip = 0;
//        for (int i = 0; i < fileLength - subSequenceLength + 1; i++) {
//            for (int j = 0; j < subSequenceLength; j++) {
//                if (shortRead[j] != refDNA[i + j])snp_flip++;
//                if (snp_flip > snp_count)break;
//            }
//            if (snp_flip <= snp_count) {
//                for (int k = 0; k < subSequenceLength; k++) {
//                    reconstructed[i + k] = shortRead[k];
//                }
//                break;
//            }
//            else {
//                snp_flip = 0;
//            }
//        }
//    }
//
//
//    reconstructed_dna << reconstructed;
//
//    RefDna.close();
//    shortReads.close();
//    reconstructed_dna.close();
//    chrono::system_clock::time_point end = chrono::system_clock::now();
//    chrono::microseconds microSec = chrono::duration_cast<chrono::microseconds>(end - start);
//    cout << "reconstruction 소요시간 : " << microSec.count() << " us\n";
//};
//
////복원된 DNA와 My Dna Sequence와의 정확도를 비교한다.
//auto rec_comparison = [](int fileLength, string recDnaPath, string myDnaPath) {
//
//    ifstream myDNA(myDnaPath);
//    ifstream recDNA(recDnaPath);
//
//    string myDnaString;
//    getline(myDNA, myDnaString);
//    string recDnaString;
//    getline(recDNA, recDnaString);
//
//    int diff = 0;
//    for (int i = 0; i < fileLength; i++) {
//        if (myDnaString[i] != recDnaString[i])diff++;
//    }
//    cout << "불일치 문자 개수 : " << diff << endl;
//    double percentage = ((double)(fileLength - diff) / (double)fileLength) * 100;
//    cout << "일치율 : " << percentage << " (%)" << endl;
//    recDNA.close();
//    myDNA.close();
//};
//
//int main() {
//    int file_length = 1000000;
//    int short_read_length = 100;
//    int num_of_snps = 11000;
//    int snp_count = 2;
//
//    cout << "테스팅 조건 Intel i7-8565U, 16GB RAM" << endl;
//    cout << "N : " << file_length << "  M : " << num_of_snps << "  L : " << short_read_length << "  D : " << snp_count << endl;
//    /*string ref_dna_file = "reference("+ to_string(file_length)+").txt";
//    string my_dna_file = "myDna("+to_string(file_length) + ", " + to_string(snp_count) + ", " + to_string(short_read_length) + ", " + to_string(num_of_snps) + ").txt";
//    string short_read_file = "shortreads("+to_string(file_length) + ", " + to_string(snp_count) + ", " + to_string(short_read_length) + ", " + to_string(num_of_snps) + ").txt";
//    string rec_dna_file = "reconstruction("+to_string(file_length) + ", " + to_string(snp_count) + ", " + to_string(short_read_length) + ", " + to_string(num_of_snps) + ")_TRI.txt";*/
//    
//    string ref_dna_file = "ref.txt";
//    string my_dna_file = "mydna.txt";
//    string short_read_file = "shortread.txt";
//    string rec_dna_file = "reconstruction(" + to_string(file_length) + ", " + to_string(snp_count) + ", " + to_string(short_read_length) + ", " + to_string(num_of_snps) + ")_TRI_Final_check.txt";
//
//    trivial_reconstruction(short_read_length, snp_count, file_length, ref_dna_file, short_read_file, rec_dna_file);
//    auto compare_DNA = [](string aDNA, string bDNA) {
//        ifstream comparisonA(aDNA);
//        ifstream comparisonB(bDNA);
//        string A; string B;
//        getline(comparisonA, A);
//        getline(comparisonB, B);
//        int count = 0;
//        for (auto i = 0; i < A.size(); i++) {
//            (A[i] != B[i]) ? ++count : NULL;
//        }
//        cout << "결과 비교 : " << 100.0 - 100 * ((double)count / (double)A.size()) << endl;
//        comparisonB.close();
//        comparisonA.close();
//    };
//
//    cout << "RefDNA-MyDNA : "; compare_DNA(ref_dna_file, my_dna_file);
//    cout << "복원 ReconstructedDNA-MyDNA : "; compare_DNA(rec_dna_file, my_dna_file);
//}