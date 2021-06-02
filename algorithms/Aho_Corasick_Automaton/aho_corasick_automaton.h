#pragma once
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

int FIND = 0;
class Component {
public:
	int length = 0;
	char twig;
	bool visited = false;
	string shortread = "";
	vector<Component*> arc;
	Component* fail_arc;
	Component() : twig('\0') { arc = {}; fail_arc = {}; }
	Component(char a) : twig(a) { arc = {}; fail_arc = {}; }
	Component* find(char a) {
		for (auto i : arc) {
			if (i->twig == a) return i;
		}
		return nullptr;
	}
}; 
vector<Component*> finalStates;
class automata {
public:

	Component* root;
	automata() {
		root = new Component('$');
	}
	void make_automaton(string short_read) {
		Component* temp = root;
		int read_length = short_read.size();
		//cout << temp->twig << endl;
		for (int i = 0; i < read_length; i++) {
			temp->length = i;
			Component* ptr = temp->find(short_read[i]);
			if (ptr == nullptr) {
				Component* newnode = new Component(short_read[i]);

				// change
				// final state info update
				// shortread : shortread string info, length : read_length (shortread size)
				if (i == read_length - 1) {
					newnode->shortread = short_read;
					newnode->length = read_length;
					//cout << finalStates.size() << endl;
					finalStates.push_back(newnode);
				}

				temp->arc.push_back(newnode);
				temp = newnode;
				//cout << temp->twig << endl;
			}
			else if (ptr != nullptr) {
				//cout << "존재함" << endl;
				temp = ptr;
			}
		}

	}
	void make_failure_route() {
		Component* temp = root;
		deque<Component*> deck = { temp };
		while (!deck.empty()) {

			Component* ptr = deck.front(); deck.pop_front(); //부모노드결정

			for (auto i : ptr->arc) { //자식노드 순회
				//#1 최초노드로부터 오토마타로 인식 실패 발견시 루트로 돌아간다.
				if (ptr == temp) {
					ptr->fail_arc = temp;
					for (auto i : ptr->arc) {
						i->fail_arc = temp;
					}
				}
				else {
					//#2 부모노드의 실패링크로 이동하고, 자식 중 접미부가 같은 문자인 것을 찾는다.
					// 못 찾을 경우 이동된 노드에서의 실패링크로 이동한다.
					Component* fail_moved = ptr->fail_arc;
					Component* needtodecide;
					do {
						needtodecide = fail_moved->find(i->twig);
						fail_moved = fail_moved->fail_arc;
					} while (fail_moved != temp && needtodecide == nullptr);

					if (needtodecide != nullptr) {
						fail_moved = needtodecide;
					}
					i->fail_arc = fail_moved;
					//cout << i->twig<<"#" << fail_moved->twig << "#"<< endl;
				}
				deck.push_back(i);
			}
		}
	}


	void print_automata(Component* temp) {

		for (int i = 0; i < temp->arc.size(); i++) {
			cout << temp->arc[i]->twig << " ";
		} cout << endl;
		for (int i = 0; i < temp->arc.size(); i++) {
			if (temp->arc[i] != nullptr)print_automata(temp->arc[i]);
		}


	}
	void print_node(Component* temp) {
		cout << "index : " << temp->length << endl;
		cout << "string : " << temp->shortread << "ShortRead Length : " << temp->shortread.length() << endl;
		cout << "Twig :" << temp->twig << endl;
		cout << "연결된 노드 :";
		for (auto i : temp->arc) {
			cout << i->twig << " ";
		}cout << endl;
		cout << "실패 아크 :";
		cout << temp->fail_arc->twig << endl;


	}

};