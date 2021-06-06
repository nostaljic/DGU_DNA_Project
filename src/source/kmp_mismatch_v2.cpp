#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<chrono>
#include<time.h>
#include<random>
#include<algorithm>
#include <random>
#include <ctime>
#include <functional>
using namespace std;
vector<char> nucleic_sequence = { 'A','T','G','C' };
vector<int> randomSeed = { 7,13,17,23 };
/*
일반적인 랜덤함수를 쓰게 될 경우 문자열(DNA)의 반복 정도가 심해집니다.
반복을 줄일 수 있는방법으로 TWIG를 int seed 값으로 모듈러 연산을 한 값이 0이 되는 순간 섞어줍니다.
그 후 seed를 randomSeed 인트 벡터의 숫자를 더해 업데이트 해줍니다.
*/
auto make_string_file = [](int fileLength, string filePath) {
    ofstream dna(filePath);
    string fileDnaString;
    random_device rd;
    int seed = 1007;
    for (int i = 0; i < fileLength; i++) {
        fileDnaString += nucleic_sequence[(rand() % 4)];
        if (i % seed == 0) {
            shuffle(nucleic_sequence.begin(), nucleic_sequence.end(), default_random_engine(rd()));
            seed += randomSeed[rand() % 4];
        }
    }
    dna << fileDnaString;
    dna.close();
};
auto endo_nuclease = [](int subSequenceLength, int shortReadsCounter, int fileLength, string filePath, string shortReadsPath) {
    vector<string> shortreads_container;
    ifstream dna(filePath);
    string fileDnaString;
    getline(dna, fileDnaString);

    ofstream shortReads(shortReadsPath);

    mt19937 engine((unsigned int)time(NULL));                    // 松本眞のアルゴリズム
    uniform_int_distribution<int> distribution(0, fileLength - subSequenceLength - 1);       // 생성 범위
    auto generator = bind(distribution, engine);
    //int maxind = fileLength - subSequenceLength;
    for (int i = 0; i < shortReadsCounter; i++) {
        int randomindex = generator();
        //int randomindex = rand()% maxind;
        shortReads << fileDnaString.substr(randomindex, subSequenceLength); shortReads << '\n';
        shortreads_container.push_back(fileDnaString.substr(randomindex, subSequenceLength));
    }

    shortReads.close();
    dna.close();
};
auto make_my_dna = [](int subSequenceLength, int snp_count, int fileLength, string refDnaPath, string myDnaPath) {

    ifstream refDNA(refDnaPath);
    string refDnaString;
    getline(refDNA, refDnaString);

    ofstream myDNA(myDnaPath);

    for (int i = 0; i < fileLength - subSequenceLength - 1; i += subSequenceLength) {
        for (int j = 0; j < snp_count - 1; j++) {
            refDnaString[i + (rand() % subSequenceLength)] = nucleic_sequence[(rand() % 4)];
        }
    }
    myDNA << refDnaString;

    refDNA.close();
    myDNA.close();
};
auto maxPrefixTable = [](string matchingString) {
    int patternLen(matchingString.length());
    vector<int> mst(patternLen, 0);   // 벡터를 만들고 0으로 초기화한다.

    for (int i(1), j(0); i < patternLen; ++i) {
        /* 1. j가 0보다 크다, 그리고 prefix와 postfix가 다르다면
        j의 뒤로 이동 및 건너뜀. */
        while (j > 0 && matchingString[i] != matchingString[j]) j = mst[j - 1];
        /* 2. prefix와 postfix가 같으면 해당 칸만 증감연산 및 한칸 뒤로 이동 */
        if (matchingString[i] == matchingString[j]) mst[i] = ++j;
    }
    return mst;   // 생성된 최대 테이블 반환한다.
};

auto KMP_matcher = [](string refFilePath, string recFilePath,string shortReadFile, int maxMismatch) {
    chrono::system_clock::time_point start = chrono::system_clock::now();
    ifstream refDna(refFilePath);
    ifstream shortRead(shortReadFile);
    ofstream reconstructedDNA(recFilePath);

    vector<int> patter;
    string ref_dna_string;
    string a_short_read;
    string my_dna_string;
    getline(refDna, ref_dna_string);

    int n = ref_dna_string.length(); 
    string reconstructed;
    reconstructed.resize(n, '\0');

    int cnt_mismatch = 0;                            // mismatch counting 변수 
    int first_mismatch_index = 0;
    int time_check = 0; 
    
    while (getline(shortRead, a_short_read)) {
        int m = a_short_read.length();
       if (time_check++ % 1000 == 0)cout << (double)time_check / 400000.0 << "\n";
        
        vector<int> mpt(maxPrefixTable(a_short_read));

        for (size_t i(0), j(0); i < n- m; ++i) {
            
            while (j > 0 && i < n - m&& ref_dna_string[i] != a_short_read[j]) {
                //cout << "hi1" << endl;
               cnt_mismatch++;  // 패턴의 j 번째에서, mismatch가 발생했음을 저장 
                if (cnt_mismatch > maxMismatch) {
                    //mismatch가 허용 개수 넘으면 
                    //j는 처음으로 발생했던 mismatch의 index - 1로 이동
                    // j = mpt[j - 1]; //mismatch 고려 안할때 이런식으로 이동
                    i = i - (j - mpt[(first_mismatch_index - 1)]) + 1;
                    j = mpt[(first_mismatch_index - 1)];

                    cnt_mismatch = 0;                  //초기화
                    first_mismatch_index = 0;
                }//TCCAGGACCGCGAGACACTTCCAACAAGGA
                 //GCCCACACGGACGCAACGCGAATAATTAGA
                ++i; ++j;
                //cout << "hi2" << endl;
            }
            if (i<n-m&&ref_dna_string[i] == a_short_read[j]) {
                //cout << "hi3" << endl;
                if (j == m - 1) {
                    //cout << i << " , " << j << " : " << a_short_read << " : " << index_vector_mismatch[0] << endl;
                    //cout << i << " , " << j << " : " << a_short_read << " : " << index_vector_mismatch[0] << endl;
                    //cout << "matched" << endl;
                    patter.push_back(i - m + 1);
                    for (int strrec(0); strrec < m; strrec++) {
                        reconstructed[i - m + strrec+1] = a_short_read[strrec];
                    }
                    j = mpt[j];
                    cnt_mismatch = 0;                  //초기화
                    first_mismatch_index = 0;     //초기화 
                    break;
                }
                else j++;
            }
        }
    }

    reconstructedDNA << reconstructed;
    refDna.close();
    shortRead.close();
    reconstructedDNA.close();
    chrono::system_clock::time_point end = chrono::system_clock::now();
    chrono::microseconds microSec = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "KMP 방식 소요시간 : " << microSec.count() << " us\n";

    return patter;
};




int main() {

    srand((unsigned)time(NULL));
    string ref_dna_file = "dna_100000000.txt";
    string my_dna_file = "snpdna_100000000.txt";
    string rec_dna_file = "reconstruced100000000DNA.txt";
    string short_read_file = "shortRead100000000.txt";
    int short_read_length = 64;
    int snp_count= 2;
    int num_of_snps = 400000;
    int file_length = 10000000;
    //100000 2000 3, 1000000 20000 147초, 71퍼센트 1000000 40000 285초, 71퍼센트
    make_string_file(file_length, ref_dna_file);
    make_my_dna(short_read_length, snp_count, file_length, ref_dna_file, my_dna_file);
    endo_nuclease(short_read_length, num_of_snps, file_length, my_dna_file, short_read_file);
    cout << "파일생성완료" << endl;
    KMP_matcher(ref_dna_file, rec_dna_file, short_read_file, 2);                               //mismatch 2개? 
    
    auto compare_DNA = [](string aDNA, string bDNA) {
        ifstream comparisonA(aDNA);
        ifstream comparisonB(bDNA);
        string A; string B;
        getline(comparisonA, A);
        getline(comparisonB, B);
        int count = 0;
        for (auto i = 0; i < A.size(); i++) {
            (A[i] != B[i]) ? ++count : NULL;
        }
        cout << "결과 비교 : " << 100.0-100*((double)count / (double)A.size()) << endl;
        comparisonB.close();
        comparisonA.close();
    };
    cout << "RefDNA-MyDNA : "; compare_DNA(ref_dna_file, my_dna_file);
    cout << "ReconstructedDNA-MyDNA : "; compare_DNA(rec_dna_file, my_dna_file);
    cout << "ReconstructedDNA-RefDNA : "; compare_DNA(rec_dna_file, ref_dna_file);

    //printPatters(filePath1, str15, resultPath1);

}
