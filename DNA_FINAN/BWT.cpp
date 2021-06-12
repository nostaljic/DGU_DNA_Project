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

#define MAX_N 1000001
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
    n = (int)strlen(str); // check num of string

    for (int i = 0; i < n; i++) // check the first char
    {
        T[i] = i;
        com_arr[i] = str[i] - 'a';
    }

    while (t <= n)
    {
        com_arr[n] = -1;
        sort(T, T + n, compare);
        com_tem_arr[T[0]] = 0; // 다음 그룹 할당 위한 번호

        for (int i = 1; i < n; i++)
        {
            if (compare(T[i - 1], T[i])) // 그룹이 다를 경우
            {
                com_tem_arr[T[i]] = com_tem_arr[T[i - 1]] + 1;
            } // 다음 그룹 번호 할당

            else // 같은 그룹인 경우
            {
                com_tem_arr[T[i]] = com_tem_arr[T[i - 1]];
            } // 해당 그룹 번호 할당
        }

        for (int i = 0; i < n; i++) // 새로운 그룹
        {
            com_arr[i] = com_tem_arr[i];
        }

        t <<= 1;
    }
}


int findSequence(int start, char bwt[], char sread[], int size, int num2, int snp_count)
{
    int copy = num2;
    int startingPoint = 0; // starting index
    int endingPoint = 0; // ending index
    int point = size - 1;

    // starting index와 ending index 지정
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
    } // checking mismatch 

    if (pointer.size() == size || bwt[start] == '$' || mis > snp_count)
    {
        return fir_vect[pri].second - 1;
    } // pointer의 사이즈가 strlen - 1와 같거나 mis가 설정한 snp_count보다 클 경우 반환
    int check1 = 0; // 초기화
    int check2 = 0; // 초기화
    check1 = bwt_vect[start].second; // 해당 문자가 동일한 문자들 중 몇번째 문자인지 마킹

    // starting point, ending point
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
    start = startingPoint + check1 - 1; // start 재정의 

    num2--;
    findSequence(start, bwt, sread, size, num2, snp_count); // 재귀
}

void find_loc(char ref[], char sread[], int size, int num2, char bwt[], int snp_count)
{
    int sta = 0;
    int check = 0;
    int point = size - 1;
    int startingPoint = 0;
    int endingPoint = 0;

    // short read의 마지막 문자에 따라 starting point, ending point 설정
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


    for (int i = startingPoint; i <= endingPoint; i++) { // starting index, ending index 조회
        if (fir_vect[i].first == sread[point]) {
            sta = i;
            pri = i;
            int num2 = size - 1;
            int mis = 0;
            check = findSequence(sta, bwt, sread, size, num2, snp_count);
            if (pointer.size() == size && mis <= snp_count) // short read가 매칭이 된 경우
            {
                int j = 0;
                for (int k = check; k < check + size; k++)
                {
                    reconstruct_BWT[k] = sread[j]; // 해당 위치의 short read를 복원
                    j++;
                }
                int sizeOfRoute = pointer.size();
                for (int k = 0; k < sizeOfRoute; k++)
                {
                    pointer.pop_back();
                } // pointer 벡터 초기화
            }
            else // 찾지 못한 경우 복원을 하지않고 벡터 초기화
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
    int short_read_length = 100;
    int snp_count = 2;
    int num_of_snps = 11000;
    int file_length = 1000000;

    cout << "테스팅 조건 Intel i7-8565U, 16GB RAM" << endl;
    cout << "N : " << file_length << "  M : " << num_of_snps << "  L : " << short_read_length << "  D : " << snp_count << endl;
   /* string ref_dna_file = "reference(" + to_string(file_length) + ").txt";
    string my_dna_file = "myDna("+to_string(file_length) + ", " + to_string(snp_count) + ", " + to_string(short_read_length) + ", " + to_string(num_of_snps) + ").txt";
    string short_read_file = "shortreads("+to_string(file_length) + ", " + to_string(snp_count) + ", " + to_string(short_read_length) + ", " + to_string(num_of_snps) + ").txt";
    string rec_dna_file = "reconstruction("+to_string(file_length) + ", " + to_string(snp_count) + ", " + to_string(short_read_length) + ", " + to_string(num_of_snps) + ")_BWT.txt";*/

    string ref_dna_file = "ref.txt";
    string my_dna_file = "mydna.txt";
    string short_read_file = "shortread.txt";
    string rec_dna_file = "reconstruction(" + to_string(file_length) + ", " + to_string(snp_count) + ", " + to_string(short_read_length) + ", " + to_string(num_of_snps) + ")_BWT_Final_check.txt";

    // BWA Algorithm
    FILE* refDNA = fopen(ref_dna_file.c_str(), "r");
    FILE* short_read_BWT = fopen(short_read_file.c_str(), "r");
    FILE* myDNA_BWT = fopen(my_dna_file.c_str(), "r");
    FILE* RE_construct_BWT = fopen(rec_dna_file.c_str(), "w");
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
    } // reconstruct 배열을 후에 맵핑되지 않은 곳을 쉽게 확인하기 위해 *로 채워준다

    strcat(ref_string, "$"); // ref dna 마지막에 $ 추가

    getT(ref_string); // suffix array 생성

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
    } // A, T, C ,G 의 갯수를 count


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
    } // BWT에 원소들을 넣어준다

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
    } // vector bwt에 배열의 원소와 동일한 문자 중 몇번째 위치에 있는지를 의미하는 정수값을 pair로 넣는다

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
    } // myDNA와 비교
    int count2 = 0;
    for (int i = 0; i < (signed)strlen(ref_string) - 1; i++)
    {
        if (reconstruct_BWT[i] != ref_string[i])
        {
            count2++;
        }
    } // refDNA와 비교

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