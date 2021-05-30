//////////////////////////////////////////////////////
////											 /////
////	 동국대학교 컴퓨터공학과				 /////
////	 수업 : 컴퓨터알고리즘과실습			 /////
////	 팀명 : 고통유전자 						 /////
////	 이름 : 박재용,윤희창,이호준,이가영		 /////
////	 <scorpion@dgu.ac.kr>					 /////
////											 /////
////	 Aho Corasick Automaton					 /////
////											 /////
//////////////////////////////////////////////////////

#include <iostream>
#include <random>
#include <ctime>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <functional>
#include <deque>
#include <map>
using namespace std;
vector<char> TWIG = { 'A','T','G','C' };
vector<string> short_reads_container;
vector<int> randomSeed = { 7,13,17,23 };
class Component {
public:
	char twig;
	vector<Component*> arc;
	vector<Component*> fail; // 실패 링크 ⭐
	Component() : twig('\0') { arc = {}; fail = {}; }
	Component(char a) : twig(a) { arc = {}; fail = {}; }
	Component* find(char a) {
		for (auto i : arc) {
			if (i->twig == a) return i;
		}
		return nullptr;
	}
};
class automata {
public:
	string maxPrefix;
	Component* root;
	automata() {
		root = new Component('$');
	}
	void make_automaton(string short_read) {
		Component* temp = root;
		int read_length = short_read.size();
		//cout << temp->twig << endl;
		for (int i = 0; i < read_length; i++) {
			Component* ptr = temp->find(short_read[i]);
			if (ptr == nullptr) {
				Component* newnode = new Component(short_read[i]);
				temp->arc.push_back(newnode);
				temp = newnode;
				//cout << temp->twig << endl;
			}
			else if (ptr != nullptr) {
				cout << "존재함" << endl;
				temp = ptr;
			}
		}

	}

	void print_automata(Component* temp) {
	
		for (int i = 0; i < temp->arc.size(); i++) {
			cout << temp->arc[i]->twig << " ";
		} cout<< endl;
		for (int i = 0; i < temp->arc.size(); i++) {
			if (temp->arc[i] != nullptr)print_automata(temp->arc[i]);
		}
		

	}

};
deque<char> aho_corasick_automata = {};

//난수 생성기 (松本眞のアルゴリズム)
auto advanced_rand = [](int range) {
	random_device rd;
	mt19937_64 engine(rd());
	uniform_int_distribution<int> distribution(0, range - 1);
	auto generator = bind(distribution, engine);
	return generator();
};
auto make_string_file = [](int fileLength, string filePath) {
	ofstream dna(filePath);
	string fileDnaString;
	random_device rd;
	int percent = 0;
	int seed = 1007;
	for (int ligation_counter = 0; ligation_counter < fileLength; ligation_counter++) {
		fileDnaString += TWIG[(rand() % 4)];
		if (ligation_counter % seed == 0) {
			shuffle(TWIG.begin(), TWIG.end(), default_random_engine(rd()));
			seed += randomSeed[rand() % 4];
		}

		
		if (ligation_counter % (fileLength / 100) == 0) {
			++percent;
			cout << percent << "% " << endl;
		}
	}
	dna << fileDnaString;
	dna.close();
};
auto make_my_dna = [](int subSequenceLength, int short_reads_count, int snp_count, int fileLength, string refDnaPath, string myDnaPath) {

	ifstream refDNA(refDnaPath);
	string refDnaString;
	getline(refDNA, refDnaString);

	ofstream myDNA(myDnaPath);

	for (int i = 0; i < fileLength - subSequenceLength - 1; i += subSequenceLength) {
		for (int j = 0; j < snp_count - 1; j++) {
			refDnaString[i + (rand() % subSequenceLength)] = TWIG[(rand() % 4)];
		}
	}
	myDNA << refDnaString;

	refDNA.close();
	myDNA.close();
};
auto endo_nuclease = [](int subSequenceLength, int shortReadsCounter, int fileLength, string filePath, string shortReadsPath) {

	ifstream dna(filePath);
	string fileDnaString;
	getline(dna, fileDnaString);

	ofstream shortReads(shortReadsPath);

	for (int i = 0; i < shortReadsCounter; i++) {
		int randomindex = advanced_rand(fileLength - subSequenceLength);
		shortReads << fileDnaString.substr(randomindex, subSequenceLength); shortReads << '\n';
	}

	shortReads.close();
	dna.close();
};

//Short reads 파일을 벡터 컨테이너에 적재
auto import_short_reads = [](string filePath) {
	ifstream short_reads_file(filePath);
	string short_read;
	while (getline(short_reads_file, short_read)) {
		short_reads_container.push_back(short_read);
	}

	for (auto i : short_reads_container) {
		cout << i << endl;
	}
};

int main() {
	
	//make_string_file(10000, "helloDna.txt");
	//make_my_dna(30,500,2,10000,"helloDna.txt", "refDna.txt");
	//endo_nuclease(30, 500, 10000, "helloDna.txt", "shortreads.txt");
	//import_short_reads("shortreads.txt");
	automata a;
	a.make_automaton("hi");
	a.make_automaton("ha");
	a.make_automaton("haha");
	a.make_automaton("hihi");
	a.make_automaton("ahihi");
	a.print_automata(a.root);
	return 0;
}