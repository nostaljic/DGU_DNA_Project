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
#define _CRT_SECURE_NO_WARNINGS
#pragma warning(disable:4996)

using namespace std;

#define MAX_N 1000000
vector<char> nucleic_sequence = { 'A','T','G','C' };
vector<int> randomSeed = { 7,13,17,23 };

// for BWT Algorithm
int pri; 
vector<char> pointer; 
vector<pair<char, int> > fir_vect; 
vector<pair<char, int> > bwt_vect; 


int t, n, com_arr[MAX_N], com_tem_arr[MAX_N], T[MAX_N];
int num_A = 0;
int num_C = 0;
int num_G = 0;
int num_T = 0;
int num_A2 = 0;
int num_C2 = 0;
int num_G2 = 0;
int num_T2 = 0;
char* ref_string = new char[MAX_N];
char reconstruct_BWT[MAX_N];
// int misMatch = 0;
char tmp[MAX_N];

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
    //vector<string> shortReads;
    ifstream dna(filePath);
    string fileDnaString;
    getline(dna, fileDnaString);

    ofstream shortReads(shortReadsPath);
    int num_of_short_reads = 0;
    int sign_flag = 1;
    int index_of_sub_dna = 0;
    //int rand_degree = 1;
    int delim = fileLength / shortReadsCounter;
    while (num_of_short_reads < shortReadsCounter - 1) {
        if (num_of_short_reads != 0)index_of_sub_dna += sign_flag * ((rand() % (subSequenceLength - delim)) + delim);
        if (index_of_sub_dna + subSequenceLength > fileLength - 1) { sign_flag = -1;  continue; }
        if (index_of_sub_dna < 0) { sign_flag = 1;  continue; }
        shortReads << fileDnaString.substr(index_of_sub_dna, subSequenceLength); shortReads << '\n';
        shortReads << fileDnaString.substr(fileLength - (subSequenceLength + index_of_sub_dna), subSequenceLength); shortReads << '\n';
        num_of_short_reads += 1;
    }
    shortReads << fileDnaString.substr(fileLength - subSequenceLength, subSequenceLength); shortReads << '\n';
    //shortreads_container.push_back(fileDnaString.substr(fileLength - subSequenceLength, subSequenceLength));

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

auto KMP_matcher = [](string refFilePath, string recFilePath, string shortReadFile, int maxMismatch) {
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
    int time_check = 0;
    int first_mismatch_index = 0;                   // 첫번째 mismatch 발생 index
    int short_read_index = 0;                       // n번째 short read에서 복원할때 사용 안된경우 나중에 뽑아내려고 쓰는 변수 
    vector<string> non_visited_shortread_vector;

    /*
    //결국 matching 못해서 복원할 때 못쓰이면
                    non_visited_shortread_vector.push_back(a_short_read);
    */
    while (getline(shortRead, a_short_read)) {

        non_visited_shortread_vector.push_back(a_short_read);

        int m = a_short_read.length();
        vector<int> mpt(maxPrefixTable(a_short_read));

        for (size_t i(0), j(0); i < n - m; ++i) {
            while (j > 0 && ref_dna_string[i] != a_short_read[j]) {
                cnt_mismatch++;
                if (cnt_mismatch == 1) first_mismatch_index = j; //첫번째 mismatch index 기억 

                if (cnt_mismatch > maxMismatch) { //mismatch가 허용 개수 넘으면 
                    //j는 처음으로 발생했던 mismatch의 index - 1로 이동
                    // j = mpt[j - 1]; //mismatch 고려 안할때 이런식으로 이동
                    //i = i - (j - mpt[(index_vector_mismatch[0] - 1)]) + 1;
                    //j = mpt[(index_vector_mismatch[0] - 1)];
                    i = i - (j - mpt[(first_mismatch_index - 1)]) + 1;
                    j = mpt[(first_mismatch_index - 1)];

                    cnt_mismatch = 0;                  //초기화
                    first_mismatch_index = 0;
                }
                ++i; ++j;

            }
            if (ref_dna_string[i] == a_short_read[j]) { //스트링 매칭 후 복원 
                if (j == m - 1) {
                    patter.push_back(i - m + 1);
                    for (int strrec(0); strrec < m; strrec++) {
                        reconstructed[i - m + strrec + 1] = a_short_read[strrec];
                    }
                    j = mpt[j];
                    cnt_mismatch = 0;                  //초기화
                    first_mismatch_index = 0;          //초기화
                    non_visited_shortread_vector.pop_back(); //복원했으면 non visited vector에서 삭제 
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

    cout << "복원할 때 안쓰인 short read의 갯수 : " << non_visited_shortread_vector.size() << endl;
    /*
    for (int i = 0; i < non_visited_shortread_vector.size(); i++) {
        cout <<non_visited_shortread_vector[i]<< endl;
    }
    */

    return non_visited_shortread_vector;
};

// BWT Algorithm

bool compare(int x, int y)
{
    if (com_arr[x] == com_arr[y])
    {
        return com_arr[x + t] < com_arr[y + t];
    }
    return com_arr[x] < com_arr[y];
}

void getT(const char* str)
{
    t = 1;
    n = (int)strlen(str);

    for (int i = 0; i < n; i++)
    {
        T[i] = i;
        com_arr[i] = str[i] - 'a';
    }

    while (t <= n)
    {
        com_arr[n] = -1;
        sort(T, T + n, compare);
        com_tem_arr[T[0]] = 0;

        for (int i = 1; i < n; i++)
        {
            if (compare(T[i - 1], T[i]))
            {
                com_tem_arr[T[i]] = com_tem_arr[T[i - 1]] + 1;
            }

            else
            {
                com_tem_arr[T[i]] = com_tem_arr[T[i - 1]];
            }
        }

        for (int i = 0; i < n; i++)
        {
            com_arr[i] = com_tem_arr[i];
        }

        t <<= 1;
    }
}


int findSequence(int start, char bwt[], char sread[], int size, int num2, int snp_count)
{
    int copy = num2;
    int startingPoint = 0;
    int endingPoint = 0;
    int point = size - 1;

    if (sread[point] == 'A')
    {
        startingPoint = 1;
        endingPoint = num_A;
    }
    else if (sread[point] == 'C')
    {
        startingPoint = num_A + 1;
        endingPoint = startingPoint + num_C - 1;
    }
    else if (sread[point] == 'G')
    {
        startingPoint = num_A + num_C + 1;
        endingPoint = startingPoint + num_G - 1;
    }
    else if (sread[point] == 'T')
    {
        startingPoint = num_A + num_C + num_G + 1;
        endingPoint = startingPoint + num_T - 1;
    }

    pointer.push_back(fir_vect[start].first);
    int mis = 0;

    for (int i = 0; i < (signed)pointer.size(); i++)
    {
        if (pointer[i] != sread[size - i - 1])
        {
            mis++;
        }
    }

    if (pointer.size() == size || bwt[start] == '$' || mis > snp_count)
    {
        return fir_vect[pri].second - 1;
    }
    int check1 = 0;
    int check2 = 0;
    check1 = bwt_vect[start].second;

    if (bwt[start] == 'A')
    {
        startingPoint = 1;
        endingPoint = num_A;
    }
    else if (bwt[start] == 'C')
    {
        startingPoint = num_A + 1;
        endingPoint = startingPoint + num_C - 1;
    }
    else if (bwt[start] == 'G')
    {
        startingPoint = num_A + num_C + 1;
        endingPoint = startingPoint + num_G - 1;
    }
    else if (bwt[start] == 'T')
    {
        startingPoint = num_A + num_C + num_G + 1;
        endingPoint = startingPoint + num_T - 1;
    }

    pri = start;
    start = startingPoint + check1 - 1;

    num2--;
    findSequence(start, bwt, sread, size, num2, snp_count);
}

// 위치 찾는 함수, 구현중
void find_loc(char ref[], char sread[], int size, int num2, char bwt[], int snp_count)
{
    int sta = 0;
    int check = 0;
    int point = size - 1;
    int startingPoint = 0;
    int endingPoint = 0;

    if (sread[point] == 'A')
    {
        startingPoint = 1;
        endingPoint = num_A;
    }
    else if (sread[point] == 'C')
    {
        startingPoint = num_A + 1;
        endingPoint = startingPoint + num_C - 1;
    }
    else if (sread[point] == 'G')
    {
        startingPoint = num_A + num_C + 1;
        endingPoint = startingPoint + num_G - 1;
    }
    else if (sread[point] == 'T')
    {
        startingPoint = num_A + num_C + num_G + 1;
        endingPoint = startingPoint + num_T - 1;
    }


    for (int i = startingPoint; i <= endingPoint; i++)
    {
        if (fir_vect[i].first == sread[point])
        {
            sta = i;
            pri = i;
            int num2 = size - 1;
            int mis = 0;
            check = findSequence(sta, bwt, sread, size, num2, snp_count);
            if (pointer.size() == size && mis <= snp_count)
            {
                int j = 0;
                for (int k = check; k < check + size; k++)
                {
                    reconstruct_BWT[k] = sread[j];
                    j++;
                }
                int sizeOfRoute = pointer.size();
                for (int k = 0; k < sizeOfRoute; k++)
                {
                    pointer.pop_back();
                }
            }
            else
            {
                int sizeOfRoute = pointer.size();
                for (int k = 0; k < sizeOfRoute; k++)
                {
                    pointer.pop_back();
                }
            }
        }
    }
}



int main() {

    srand((unsigned)time(NULL));
    string ref_dna_file = "dna_100000000.txt";
    string my_dna_file = "snpdna_100000000.txt";
    string rec_dna_file = "reconstruced100000000DNA.txt";
    string short_read_file = "shortRead100000000.txt";

    int short_read_length = 64;
    int snp_count = 2;
    int num_of_snps = 100000;
    int file_length = 2000000;

    make_string_file(file_length, ref_dna_file);
    make_my_dna(short_read_length, snp_count, file_length, ref_dna_file, my_dna_file);
    endo_nuclease(short_read_length, num_of_snps, file_length, my_dna_file, short_read_file);
    cout << "파일생성완료" << endl << endl;

    ofstream unvisitedShortReads("unvisited_file.txt");
    for (auto i : KMP_matcher(ref_dna_file, rec_dna_file, short_read_file, 2)) {
        unvisitedShortReads << i << '\n';
    } // mismatch 2
    // KMP Algorithm
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
        cout << "결과 비교 : " << 100.0 - 100 * ((double)count / (double)A.size()) << "%" << endl;
        comparisonB.close();
        comparisonA.close();
    };

    cout << "RefDNA-MyDNA : "; compare_DNA(ref_dna_file, my_dna_file);
    cout << "ReconstructedDNA-MyDNA : "; compare_DNA(rec_dna_file, my_dna_file);
    cout << "ReconstructedDNA-RefDNA : "; compare_DNA(rec_dna_file, ref_dna_file);
    //printPatters(filePath1, str15, resultPath1);

    cout << endl;

    // BWA Algorithm
    FILE* refDNA = fopen("dna_100000000.txt", "r"); 
    FILE* short_read_BWT = fopen("shortRead100000000.txt", "r"); 
    FILE* myDNA_BWT = fopen("snpdna_100000000.txt", "r"); 
    FILE* RE_construct_BWT = fopen("reconstruct_BWT.txt", "w"); 
    char* BWT = new char[MAX_N];
    char* shortread_bwt = new char[MAX_N];
    char* myDNA = new char[MAX_N];
    int fullsize = 0;
    char ch;
    while ((ch = fgetc(refDNA)) != EOF)  
    {
        ref_string[fullsize] = ch;
        fullsize++;
    }
    ref_string[fullsize] = '\0';
    fclose(refDNA); 

    for (int k = 0; k < fullsize; k++)
    {
        reconstruct_BWT[k] = '*';
    }
    
    strcat(ref_string, "$"); 
    
    getT(ref_string);

    for (int i = 0; i < (signed)strlen(ref_string); i++) 
    {
        fir_vect.push_back(make_pair(ref_string[T[i]], T[i])); 
        if (ref_string[i] == 'A') 
        {
            num_A++;
        }
        if (ref_string[i] == 'C')  
        {
            num_C++;
        }
        else if (ref_string[i] == 'G')  
        {
            num_G++;
        }
        else if (ref_string[i] == 'T')  
        {
            num_T++;
        }
    }


    for (int i = 0; i < n; i++)
    {
        if (T[i] == 0)
        {
            BWT[i] = ref_string[strlen(ref_string) - 1];
        }
        else
        {
            BWT[i] = ref_string[T[i] - 1];
        }
    }

    for (int i = 0; i < n; i++) 
    {
        if (BWT[i] == 'A')
        {
            num_A2++;
            bwt_vect.push_back(make_pair(BWT[i], num_A2));
        }
        else if (BWT[i] == 'C')
        {
            num_C2++;
            bwt_vect.push_back(make_pair(BWT[i], num_C2));
        }
        else if (BWT[i] == 'G')
        {
            num_G2++;
            bwt_vect.push_back(make_pair(BWT[i], num_G2));
        }
        else if (BWT[i] == 'T')
        {
            num_T2++;
            bwt_vect.push_back(make_pair(BWT[i], num_T2));
        }
        else if (BWT[i] == '$')
        {
            bwt_vect.push_back(make_pair(BWT[i], 1));
        }
    }

    int j = 0;
    int short_length_BWT = short_read_length + 1;

    chrono::system_clock::time_point start = chrono::system_clock::now();
    while (!feof(short_read_BWT)) 
    {
        fgets(shortread_bwt, short_length_BWT + 1, short_read_BWT);
        find_loc(ref_string, shortread_bwt, short_length_BWT - 1, short_length_BWT - 2, BWT, snp_count);
        j++;
        continue;
    }
    fclose(short_read_BWT); 
    int k = 0;

    while ((ch = fgetc(myDNA_BWT)) != EOF)  
    {
        myDNA[k] = ch; 
        k++;
    }
    myDNA[k] = '\0';
    fclose(myDNA_BWT); 

    for (int i = 0; i < (signed)strlen(ref_string) - 1; i++) 
    {
        if (reconstruct_BWT[i] - '0' != 0)
        {
            fprintf(RE_construct_BWT, "%c", reconstruct_BWT[i]);
        }
    }
    fclose(RE_construct_BWT);

    int count1 = 0;
    for (int i = 0; i < (signed)strlen(ref_string) - 1; i++) 
    {
        if (reconstruct_BWT[i] != myDNA[i])
        {
            count1++;
        }
    }
    int count2 = 0;
    for (int i = 0; i < (signed)strlen(ref_string) - 1; i++) 
    {
        if (reconstruct_BWT[i] != ref_string[i])
        {
            count2++;
        }
    }
    cout << endl << endl << "BWT Algorithm" << endl;
    long double num = (((long double)strlen(ref_string) - 1) - count1) * 100.0;
    double percent = (num / (strlen(ref_string) - 1));
    cout << "ReconstructedDNA-MyDNA : 결과 비교 : " << percent << "% " << endl;
    long double num2 = (((long double)strlen(ref_string) - 1) - count2) * 100.0;
    double percent2 = (num2 / (strlen(ref_string) - 1));
    cout << "ReconstructedDNA-refDNA : 결과 비교 : " << percent2 << "% " << endl << endl;
    chrono::system_clock::time_point end = chrono::system_clock::now();
    chrono::microseconds microSec2 = chrono::duration_cast<chrono::microseconds>(end - start);
    cout << "BWT 방식 소요시간 : " << microSec2.count() << " us\n";
    delete[] ref_string; 
    delete[] BWT;     
    delete[] myDNA;   
    delete[] shortread_bwt;  
   

}