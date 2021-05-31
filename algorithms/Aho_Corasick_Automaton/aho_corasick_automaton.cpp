﻿//////////////////////////////////////////////////////
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
#include "aho_corasick_automaton.h"
auto automaton_dna_reconstruction = [](int subSequenceLength, int snp_count, int fileLength, string refFilePath, string reconstructionFilePath, automata aca) {

	ofstream reconstructed_dna(reconstructionFilePath);
	string reconstructed; reconstructed.resize(fileLength + 1);
	//cout << reconstructed << endl;
	ifstream RefDna(refFilePath);
	string refDNA;
	getline(RefDna, refDNA);


	for (int i = 0; i < fileLength - subSequenceLength; i += subSequenceLength) {
		Component* temp = aca.root;
		Component* ptr;
		string tempstr; tempstr.resize(30);
		int flag = 0;
		int ref = 0;
		for (int j = 0; j < subSequenceLength; j++) {

			do {
				ptr = temp->find(refDNA[i + j]);
				if (ptr != nullptr ) {
					tempstr[j] = ptr->twig;
					temp = ptr;
					flag++;
					break;
				}
				else {
					temp = temp->fail_arc;

				}

			} while (ptr == nullptr);
			cout << j << " # " << tempstr << endl;
			ref++;
		}
		if (flag == 30) {
		}
		cout << ref << endl;

	}
	//subSequenceLength==automata의 arc length
};



// DFS 함수
// parameter info
// current : 현재 확인할 위치
// flag : dfs 성공 여부
void inter_DFS(Component* ptr, bool visited[], int start, int subSequenceLength, string& ref, int current, int mis_match, int snip_count, bool& flag) {

	// find shortread AND update Reference
	if (ptr->length == subSequenceLength && flag == false) {
		flag = true;
		string shortread = ptr->shortread;
		int index = current - subSequenceLength;
		for (int i = 0; i < subSequenceLength; i++) {
			cout << "==============dfs change : " << ref[index + i] << "  " << shortread[i] << endl;
			ref[index + i] = shortread[i];
		}
		// 6.01 해당 shortread 방문
		ptr->visited = true;
		FIND++;
		cout << "\n";
		return;
	}

	visited[start] = true;
	int save = mis_match;

	// 다음 노드 확인
	for (int i = 0; i < ptr->arc.size(); i++) {
		Component* temp = ptr->arc[i];
		//cout << "check in dfs : " << ref[current] << " " << temp->twig <<" "<< mis_match<< endl;

		// 해당 노드 값이 reference와 다를 때 
		if (ref[current] != temp->twig) {
			mis_match++;
			if (mis_match > snip_count) {
				// 허용 mis_match를 초과했을 때 더이상 탐색 x
				mis_match = save;
				continue;
			}
		}

		// #case 1 : not find && mis_match < snip_cout
		// #case 2 : find
		if (visited[temp->length] == false) {
			//cout << "check in dfs : " << ref[current] << " " << temp->twig <<" "<< mis_match<< endl;
			inter_DFS(temp, visited, temp->length, subSequenceLength, ref, current + 1, mis_match, snip_count, flag);
		}
		mis_match = save; // mismatch정보 초기화
	}
};

// MYDNA reconstruct 함수
// 목표 : reference를 이용해 shortread 오토마타 상에서 자유롭게 돌아다니면서 shortread를 찾아내고 그 결과를 reference에 반영하여 reconstruct
// Case #1. no missmatch & find shortread     => reference에 반영
// Case #2. mis_match 존재 => 1. 허용 mis_match내에서 자식 트라이를 DFS를 통해 탐색하여 shortread찾기 2. 만약 dfs로 못찾을 경우 faillink이용해 jump
// Case #3. mis_match 허용치 초과 => faillink이용해 새로운 탐색시작
// Case #4. 만약 root로 jump된다면 while(1)문 빠져나오고 새로운 i로 탐색 시작
auto automaton_dna_reconstruction_v1 = [](int subSequenceLength, int snp_count, int fileLength, string refFilePath, string reconstructionFilePath, automata aca) {

	ofstream reconstructed_dna(reconstructionFilePath);
	string reconstructed; reconstructed.resize(fileLength + 1);
	//cout << reconstructed << endl;
	ifstream RefDna(refFilePath);
	string reconstructDNA;
	getline(RefDna, reconstructDNA);

	int mis_match = 0; // mis_match count
	int position = 0; // reference에서 현재 확인하고 있는 DNA sequence 위치
					// ex) i + position => i부터 시작하는 reference의 현재 보는 char 위치(=position)

	for (int i = 0; i < fileLength - subSequenceLength; i++) {
		Component* temp = aca.root;
		Component* ptr = nullptr;
		string find_shortread;

		// 새로운 탐색 시작 초기화
		mis_match = 0;
		position = 0;

		// aho_corasick 오토마타 상에서 shortread 찾는 while문
		while (1) {

			if (ptr != nullptr && ptr->twig == '$') {
				// if ptr == root 다음 reference DNA index부터 시작
				// root를 만나면 break;
				break;
			}

			// find shortRead 
			if (ptr != nullptr && ptr->length == subSequenceLength) {
				find_shortread = ptr->shortread;

				// ShortRead ReferenceDNA에 update
				for (int a = 0; a < subSequenceLength; a++) {
					//cout << "test : " <<i <<" " << reconstructDNA[i + a] << " " << find_shortread[a] <<endl;
					reconstructDNA[i + a] = find_shortread[a];
				}
				//	cout <<"mis_match : "<<mis_match<<" "<<position<< endl<<endl;

				// 6. 01 노드 방문
				ptr->visited = true;
				FIND++;
				i += (subSequenceLength - 1); // break 이후 while문 나간 후 for문의 i++ 이후 i < fileLength - subSequenceLength test 하기 때문
				break;
			}


			// temp->find()를 통해 reference의 해당 노드를 오토마타 상에서 찾기
			if (i + position > fileLength - subSequenceLength) break; 			// string index overflow 차단
			ptr = temp->find(reconstructDNA[i + position]);


			// find의 해당 문자가 존재하는 경우
			if (ptr != nullptr) {
				temp = ptr;
				position++;
			}
			else if (ptr == nullptr && mis_match < snp_count) {
				// find결과 해당 문자 못찾았지만 mis_match < snp_count라 스닙일 가능성 존재하여 DFS로 자식 전체 탐색 진행

				bool* visited = new bool[subSequenceLength + 1];
				memset(visited, false, sizeof(visited));

				bool flag = false; // flag : dfs 찾으면 끝내도록

				// DFS
				inter_DFS(temp, visited, temp->length, subSequenceLength, reconstructDNA, i + position, mis_match, snp_count, flag);

				// DFS 결과 찾으면 flag == true
				if (flag == false) {
					// goto faillink ( faillink에서 새로운 탐색 )
					temp = temp->fail_arc;
					if (temp == nullptr) break;
					mis_match = 0;	// 초기화 
					i = i + position - temp->length; // i의 위치 변경 mismatch position과 fail link를 통해 위치 연산
					position = temp->length; // position을 fail link의 length로 변경해 바로 다음 노드 부터 탐색할 수 있도록 
				}
				else {
					// 내부 dfs shortread 탐색 성공
					// 다음 subSequenceLength부터 새로운 탐색 시작하기 위해 i 증가
					i += (subSequenceLength - 1);
				}
			}
			else {
				// jump to fail_arc;
				temp = temp->fail_arc;
				if (temp == nullptr) break;
				mis_match = 0;
				i = i + position - temp->length;
				position = temp->length;
			}
		}
	}
	// reconstruct 끝
	reconstructed_dna << reconstructDNA;
	return reconstructDNA;
};



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

	/*for (auto i : short_reads_container) {
		cout << i << endl;
	}*/
};

auto reconstruct_precision = [](string reconstruct, string mydna) {
	int miss_count = 0;
	int dna_size = reconstruct.size();
	for (int i = 0; i < dna_size; i++) {
		if (reconstruct[i] != mydna[i]) {
			miss_count++;
		}
	}
	cout << "MYDNA 복원결과 : " << double(dna_size - miss_count) * 100 / dna_size << "%입니다" << endl;
};


auto loadMyDna = [](string filename) {
	string mydna;
	ifstream loadDNA;
	loadDNA.open(filename);
	loadDNA >> mydna;
	return mydna;
};

int main() {
	// 30 50 2 1000
	int fileSize = 100000;
	int subSequenceLength = 30;
	int subSequenceCount = 4000;
	int snpCount = 2;
	make_string_file(fileSize, "refDna.txt");
	make_my_dna(subSequenceLength, subSequenceCount, 2, fileSize, "refDna.txt", "myDna.txt");
	endo_nuclease(subSequenceLength, subSequenceCount, fileSize, "myDna.txt", "shortreads.txt");
	import_short_reads("shortreads.txt");

	automata ACA;
	string reconstruct;
	string mydna;
	string reference;
	for (auto i : short_reads_container) {
		ACA.make_automaton(i);
	}
	ACA.make_failure_route();

	reference = loadMyDna("refDna.txt");

	ACA.print_node(ACA.root->arc[0]->arc[0]);//0 1 2 
	cout << "Automata 구성 완료" << endl;
	//automaton_dna_reconstruction(30, 2, 1000, "refDna.txt", "reconstruction.txt", ACA);
	//reconstruct = automaton_dna_reconstruction_v1(30, 3, 1000000, "helloDna.txt", "reconstruction.txt", ACA);
	reconstruct = automaton_dna_reconstruction_v1(subSequenceLength, snpCount, fileSize, "refDna.txt", "reconstruction.txt", ACA);
	mydna = loadMyDna("myDna.txt");
	cout << "My DNA vs Reference" << endl;
	reconstruct_precision(reference, mydna);

	cout << "My DNA vs reconstruct" << endl;
	reconstruct_precision(reconstruct, mydna);

	int check = 0; int visited = 0;

	for (auto& a : finalStates) {
		if (a->length != 30) {
			cout << "ERROR발생 length : " << a->length << endl;
		}
		if (a->visited) visited++;
		check++;
	}
	cout << "확인된 shortRead 개수 : " << check << endl;
	cout << "총 방문한 shortRead 개수 : " << visited << endl;

	int count = 0;
	sort(short_reads_container.begin(), short_reads_container.end());
	for (int i = 0; i < short_reads_container.size() - 1; i++) {
		if (short_reads_container[i] != short_reads_container[i + 1]) {
			count++;
		}
	}
	cout << "중복제외 shortRead 개수 : " << count << endl;
	return 0;
}