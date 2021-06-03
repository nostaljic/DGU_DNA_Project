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
char nucleic_sequence[] = { 'A','T','G','C' };

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
auto maxPrefixTable = [](string matchingString) {
	int patternLen(matchingString.length());
	vector<int> mst(patternLen, 0);	// 벡터를 만들고 0으로 초기화한다.

	for (int i(1), j(0); i < patternLen; ++i) {
		/* 1. j가 0보다 크다, 그리고 prefix와 postfix가 다르다면
		j의 뒤로 이동 및 건너뜀. */
		while (j > 0 && matchingString[i] != matchingString[j]) j = mst[j - 1];
		/* 2. prefix와 postfix가 같으면 해당 칸만 증감연산 및 한칸 뒤로 이동 */
		if (matchingString[i] == matchingString[j]) mst[i] = ++j;
	}
	return mst;	// 생성된 최대 테이블 반환한다.
};
auto KMP_matcher = [](string refFilePath, string shortReadsPath) {
	chrono::system_clock::time_point start = chrono::system_clock::now();
	ifstream refDna(refFilePath);
	ifstream shortRead("shortRead.txt");
	ofstream reconstructedDNA("mismatch0.txt");
	vector<int> patter;
	string ref_dna_string;
	string a_short_read;
	getline(refDna, ref_dna_string);

	int n = ref_dna_string.length();
	string reconstructed; reconstructed.resize(n, '\0');

	while (getline(shortRead, a_short_read)) {
		int m = a_short_read.length();
		vector<int> mpt(maxPrefixTable(a_short_read));

		for (size_t i(0), j(0); i < n; ++i) {
			while (j > 0 && ref_dna_string[i] != a_short_read[j]) j = mpt[j - 1]; if (ref_dna_string[i] == a_short_read[j]) {
				if (j == m - 1) {
					patter.push_back(i - m + 1);
					for (int strrec(0); strrec < m; strrec++) {
						reconstructed[i - m + strrec] = a_short_read[strrec];
					}
					j = mpt[j];
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
	cout << patter[0] << endl;
	return patter;
};
auto make_string_file = [](int fileLength, string filePath) {
	ofstream dna(filePath);
	string fileDnaString;
	for (int i = 0; i < fileLength; i++) {
		fileDnaString+=nucleic_sequence[(rand()%4)];
	}
	dna << fileDnaString;
	dna.close();
};


int main() {


	srand((unsigned)time(NULL));
	string filePath1 = "dna_10000.txt";

	
	/*******************파일생성*******************/
	//make_string_file(10000, filePath1);
	//make_string_file(100000, filePath2);
	//make_string_file(1000000, filePath3);
	//make_string_file(10000000, filePath4);
	/**********************************************/

	string str15 = "CCTGGTATGCGATC";



	string resultPath1 = "result_10000_15.txt";
	endo_nuclease(30, 5000, 10000, filePath1, "shortRead.txt");
	KMP_matcher(filePath1, str15);
	//printPatters(filePath1, str15, resultPath1);


}
