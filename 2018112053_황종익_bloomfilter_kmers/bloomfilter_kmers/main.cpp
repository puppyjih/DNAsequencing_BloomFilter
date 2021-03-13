// 2018112053 황종익

#include <iostream>
#include <vector>
#include <math.h>
#include <unordered_set>
#include <fstream>
#include <string>
#include <list>
#include <map>
#include <array>
#include <sstream>
#include <time.h>
#include "MurmurHash3.h"
#define shortReadLen 100
#define maxBreadth 90000
#define maxDepth 90000
using namespace std;

string maxLenDNA;

namespace Generator {
	const vector<string> BASES = { "A", "T", "C", "G" };

	vector<string> genExtensions(string kmer)
	{
		vector<string> E;

		unsigned long kmerSize = kmer.size();
		for (string b : BASES)
		{
			E.push_back(b + kmer.substr(0, kmerSize - 1));
			E.push_back(kmer.substr(1) + b);
		}
		
		return E;
	}

	vector<string> genLeftExtensions(string kmer)
	{
		vector<string> E;

		unsigned long kmerSize = kmer.size();
		for (string b : BASES)
		{
			E.push_back(b + kmer.substr(0, kmerSize - 1));
		}

		return E;
	}

	vector<string> genRightExtensions(string kmer)
	{
		vector<string> E;

		unsigned long kmerSize = kmer.size();
		for (string b : BASES)
		{
			E.push_back(kmer.substr(1) + b);
		}

		return E;
	}

	string extLastKmerInSeq(string sequence, int k)
	{
		return sequence.substr(sequence.size() - k);
	}

	string extFirstKmerInSeq(string sequence, int k)
	{
		return sequence.substr(0, k);
	}

}

namespace measures
{
	float measure_n50(vector<int> lengths)
	{
		map<int, int> freq;
		
		for (int i : lengths)
			freq[i] = freq[i] + 1;

		vector<int> tmpLengths;

		for (const auto& p : freq)
		{
			int key = p.first;
			int value = p.second;

			for (int i = 0; i < value * key; i++)
				tmpLengths.push_back(key);
		}

		size_t size = tmpLengths.size();

		if (size % 2 == 0) 
			return (tmpLengths[size / 2 - 1] + tmpLengths[size / 2]) / 2.;
		else
			return tmpLengths[size / 2];
	}

	vector<int> read_seq_counts(string inputpaths)
	{
		ifstream in(inputpaths);

		vector<int> counts;
		string delimiter = "길이: ";

		string line;
		while (getline(in, line))
		{
			if (line[0] >= '0' && line[0] <= '9')
			{
				string token = line.substr(line.find(delimiter) + delimiter.size());
				//cout << "token : " << token << "\n";
				counts.push_back(stoi(token));
			}
		}

		in.close();
		return counts;
	}
}

class BloomFilter
{
	int numHash;
	vector<bool> bits;

public:
	BloomFilter(int size, int numHash) : bits(size), numHash(numHash) {}

	BloomFilter(unsigned int mer_counts, int k)
	{
		int filterSize = calc_filter_size(mer_counts, k);
		int numHashFunctions = calc_numOf_hashfunctions(filterSize, mer_counts);
		cout << "filterSize : " << filterSize << " numHashFunctions : " << numHashFunctions << "\n";
		/*string s;
		cin >> s;*/
		// 필터 사이즈 조정
		while (filterSize % numHashFunctions != 0)
			filterSize++;

		new (this) BloomFilter(filterSize, numHashFunctions);
	}

	void add(const string s)
	{
		array<uint64_t, 2> hashes = hash(&s);
		//cout << "s : " << s << "\n";
		for (int n = 0; n < numHash; n++)
		{
			bits[mixHash(n, hashes[0], hashes[1], this->bits.size())] = true;
		}
		//cout << "\n\n";
	}

	bool contains(const string s)
	{
		array<uint64_t, 2> hashes = hash(&s);
		//cout << "contains s : " << s << "\n";
		for (int n = 0; n < this->numHash; n++)
		{
			if (bits[mixHash(n, hashes[0], hashes[1], this->bits.size())] == false)
			{
				return false;
			}
		}
		//cout << "\n";
		return true;
	}

	unsigned long byteSize()
	{
		return bits.size();
	}

private:
	array<uint64_t, 2> hash(const string* s)
	{
		array<uint64_t, 2> hashes;
		uint32_t seed = 420;
		MurmurHash3_x64_128(s->c_str(), (int)s->size(), seed, hashes.data());
		return hashes;
	}

	uint64_t mixHash(uint8_t n, uint64_t hashA, uint64_t hashB, uint64_t filterSize)
	{
		uint64_t h = (hashA + n * hashB) % filterSize;
		return h;
	}

	int calc_filter_size(int kmerSize, int k)
	{
		return (int)(1.44 * kmerSize * log2(1 / (2.08 / (16 * k))));
	}
	int calc_numOf_hashfunctions(int filterSize, float kmerSize)
	{
		return (int)filterSize / kmerSize * log(2);
	}
};

class DeBruijnGraph
{
public:
	DeBruijnGraph(string inputPath, unsigned int mer_counts, int k)
		: k(k), bloomFilter(mer_counts, k)
	{
		initBloomFilter(inputPath);
		findCriticalFP(inputPath);
	}
	
	void traversal(string inputPath, string outputPath, bool ignoreWords)
	{
		cout << "그래프 순회 중...\n";

		unordered_set<string> initKmers;
		ifstream in(inputPath);

		string line;
		//string delimiter = "번째 contig, 길이: ";
		while (getline(in, line))
		{
			if (ignoreWords)
				getline(in, line);

			//vector<string> leftExtensions = Generator::genLeftExtensions(Generator::extFirstKmerInSeq(line, k));
			vector<string> leftExtensions = Generator::genLeftExtensions(line);

			int index = 0;
			for (string e : leftExtensions)
			{
				if (isBelongsToDebruijnGraph(e))
					index++;
			}
			if (index == 0) 
			{
				initKmers.insert(line);
			}
		}

		in.close();

		unsigned long initKmersSize = initKmers.size();
		cout << "initKmersSize : " << initKmersSize << "\n";
		int kmerIndex = 1;

		unordered_set<string> contigs; 
		unordered_set<string> marked; 
		for (string start : initKmers)
		{
			cout << "Starting kmer: " << start << ", 진행도 : " << kmerIndex++ << "/" << initKmersSize << "\n";
			list<string> paths;
			paths.push_back(start);

			int depth = 0; 
			while (paths.size() > 0)
			{
				if (depth >= maxDepth)
				{
					contigs.insert(paths.back());
					break;
				}

				list<string> pathsBeAdd;
				list<string> pathsBeRemove;

				int breadth = paths.size();
				if (breadth == 1)
				{
					for (string& path : paths)
					{
						string lastKmer = Generator::extLastKmerInSeq(path, k);
						/*cout << "lastKmer : " << lastKmer << " size : " << lastKmer.size() << "\n";*/
						marked.insert(lastKmer);

						vector<string> ext = Generator::genRightExtensions(lastKmer);
						vector<string> validExt;
						for (string e : ext)
						{
							if (isBelongsToDebruijnGraph(e) && marked.count(e) == 0)
							{
								validExt.push_back(e);
							}
						}

						unsigned long validExtCnt = validExt.size();
						if (validExtCnt == 0)
						{
							contigs.insert(path);
							pathsBeRemove.push_back(path);
						}
						else if (validExtCnt == 1)
						{
							char& lastChar = validExt[0].back();
							path.push_back(lastChar);
						}
						else
						{
							contigs.insert(path);
							pathsBeRemove.push_back(path);
							for (string validE : validExt)
							{
								//cout << "pathsBeAdd from : " << path << " to : " << validE << "\n";
								pathsBeAdd.push_back(validE);
							}
							depth++;
						}
					}
				}
				else
				{
					depth++;
					for (string& path : paths)
					{
						string lastKmer = Generator::extLastKmerInSeq(path, k);
						//cout << "lastKmer : " << lastKmer << "\n";
						marked.insert(lastKmer);

						vector<string> ext = Generator::genRightExtensions(lastKmer);
						vector<string> validExt;
						for (string e : ext)
						{
							if (isBelongsToDebruijnGraph(e) && marked.count(e) == 0)
								validExt.push_back(e);
						}

						string pathClone(path);
						unsigned long validExtCnt = validExt.size();
						if (validExtCnt == 0)
						{
							pathsBeRemove.push_back(path);
						}
						else if (validExtCnt == 1)
						{
							char& lastChar = validExt[0].back();
							path.push_back(lastChar);
						}
						else
						{
							for (string validE : validExt)
							{
								char& lastChar = validE.back();
								string newPath(pathClone);
								newPath.push_back(lastChar);
								pathsBeAdd.push_back(newPath);
							}
						}
					}
				}

				for (string p : pathsBeRemove)
					paths.remove(p);

				if (paths.size() == 1)
					depth = 0;

				for (string p : pathsBeAdd)
				{
					if (paths.size() >= maxBreadth)
						break;
					paths.push_back(p);
				}
			}
		}

		ofstream output;
		output.open(outputPath);
		int outputIndex = 0;
		for (string contig : contigs)
		{
			/*if (contig.length() < 2 * shortReadLen + 1) continue;*/

			output << outputIndex++ << "번째 contig, 길이: " << contig.size() << "\n";
			if (contig.size() > maxLenDNA.size())
				maxLenDNA = contig;
			output << contig << "\n";
		}
		
		output.close();
	}

	unsigned long byteSize()
	{
		unsigned long size = 0;
		size += bloomFilter.byteSize();
		size += criticalFP.size() * k;
		return size;
	}

private:
	int k;
	BloomFilter bloomFilter;
	unordered_set<string> criticalFP; // critical false positive(거짓 긍정)

	void initBloomFilter(string inputPath)
	{
		cout << "블룸 필터 생성중...\n";
		unordered_set<string> kmers = loadKmers(inputPath);

		for (string kmer : kmers)
			bloomFilter.add(kmer);
	}

	void findCriticalFP(string inputPath)
	{
		cout << "거짓 긍정들을 찾아내고 있습니다...\n";

		unordered_set<string> S = loadKmers(inputPath);
		unordered_set<string> P = findP(S);

		for (string p : P)
		{
			if (S.find(p) == S.end())
				criticalFP.insert(p);
		}
	}

	unordered_set<string> findP(unordered_set<string>& S)
	{
		unordered_set<string> P;

		for (string s : S)
		{
			vector<string> E = Generator::genExtensions(s);
			for (string e : E)
			{
				if (bloomFilter.contains(e))
					P.insert(e);
			}
		}

		return P;
	}

	bool isBelongsToDebruijnGraph(string kmer)
	{
		//cout << kmer << " bloomFilter.contains(kmer) : " << bloomFilter.contains(kmer) << "\n";
		return bloomFilter.contains(kmer) && criticalFP.count(kmer) == 0;
	}

	unordered_set<string> loadKmers(string path)
	{
		ifstream in(path);

		unordered_set<string> S;
		string line;
		while (getline(in, line))
		{
			S.insert(line);
		}
		in.close();

		return S;
	}
};

vector<string> read_mers(string inputPath, bool ignoreWords)
{
	ifstream in(inputPath);
	vector<string> mers;
	//string delimiter = "번째 contig, 길이: ";

	string line;
	while (getline(in, line))
	{
		if(ignoreWords)
			getline(in, line);
		mers.push_back(line);
	}

	in.close();
	return mers;
}

void make_mers(string inputPath, string outputPath, int k, bool ignoreWords)
{
	ifstream in(inputPath);
	ofstream out(outputPath);
	
	string line;
	string delimiter = "번째 contig, 길이: ";
	while (getline(in, line))
	{
		int startIdx = 0;
		int lengthOfline = line.size();
		if (ignoreWords == false)
		{
			while (startIdx + k <= lengthOfline)
			{
				out << line.substr(startIdx++, k) << "\n";
			}
		}
		else
		{
			getline(in, line);
			if (k <= lengthOfline)
			{
				out << line.substr(0, k) << "\n";
				out << line.substr(line.size() - k) << "\n";
			}
		}
	}

	in.close();
	out.close();
}

int main()
{
	ios_base::sync_with_stdio(0);
	cin.tie(0);
	cout.tie(0);

	clock_t start = clock();
	int k; // k-mer의 개수
	
	string mainPath = "./";
	string inputPath = mainPath + "shortread_1000.txt";
	string myDNAPath = mainPath + "Mydna_1000.txt";
	string kmerPath = mainPath + "kmersRead.txt";
	string outputPath = mainPath + "repair_dna_1000.txt";
	bool calcN50 = true; // n50이 나오지 않는 경우에는 false를 사용하는 게 좋음.
	
	string tmp;
	ifstream in(inputPath);
	getline(in, tmp);
	if (tmp.size() > 50)
		k = 24;
	else if (tmp.size() > 25)
		k = 15;
	else if (tmp.size() > 10)
		k = 7;
	else if (tmp.size() > 0)
		k = 3;

	make_mers(inputPath, kmerPath, k, false);
	vector<string> kmers = read_mers(kmerPath, false);
	cout << "kmers의 수: " << kmers.size() << "\n";
	DeBruijnGraph graph = DeBruijnGraph(kmerPath, kmers.size(), k);
	graph.traversal(kmerPath, outputPath, false);


	clock_t end = clock();
	double time = double((double)end - start) / CLOCKS_PER_SEC;
	cout << "소모 시간 : " << time << "seconds\n";
	cout << "De Bruijn 그래프 크기(바이트) : " << graph.byteSize() << "\n";

	string myDNA;
	string restoredDNA = maxLenDNA;
	unsigned long long cnt = 0, highestCnt = 0;
	float n50_measure, accuracy = 0;
	ifstream in1(myDNAPath);
	getline(in1, myDNA);

	cout << "\n\n";
	cout << "      원본 myDNA의 길이      :" <<  myDNA.size() << "\n";
	cout << "복구된 가장 긴 contig의 길이 :" << restoredDNA.size() << "\n";

	cout << "\nString Matching 중...\n";
	for (int i = 0; i < myDNA.size() - restoredDNA.size(); i++)
	{
		cnt = 0;
		for (int j = 0; j < restoredDNA.size(); j++)
		{
			if (restoredDNA[j] == myDNA[i + j])
				cnt++;
		}
		if (highestCnt < cnt)
			highestCnt = cnt;
	}
	cout << "String Matching 완료!\n";

	accuracy = double(highestCnt / (double)myDNA.size()) * 100.0f;
	//highestCnt = myDNA.find(restoredDNA);
	//if (highestCnt < myDNA.size())
		

	cout << "\n다음에 나오는 아래의 정확도는 원본 길이의 1/4배정도가 복구되었을 때만 의미가 있는 값입니다.\n그보다 짧다면 더 아래의 N50값을 확인해주세요.\nDNA 일치율 : " << accuracy << "%\n";
	if (calcN50 == true)
	{
		vector<int> lengths = measures::read_seq_counts(outputPath);
		n50_measure = measures::measure_n50(lengths);
		cout << "N50 : " << n50_measure << "\n";
	}

	if (restoredDNA.size() >= myDNA.size() / 4)
	{
		cout << "(유의미한 경우에 속한다.) DNA 일치율 : " << accuracy << "%\n";
	}
	else
	{
		cout << "(유의미한 경우에 속한다.) N50의  값  : " << n50_measure << "\n";
	}

}