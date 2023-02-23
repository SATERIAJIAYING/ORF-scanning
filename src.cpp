#include <iostream>
#include <string>
#include <algorithm>
#include <assert.h>
#include <list>
#include <math.h>
#include "src.h"


int main()
{
	const std::string s1 = "atcgaagtatcggagagaatgaaaaaaaagtaaaaataacgtataagg";
	OrfScanning(std::cout, s1, 10);
	OrfScanning(std::cout, s1, 20);
	const std::string s2 = "GCAGGTACATGGCCAAGTACTAAAACG";
	OrfScanning(std::cout, s2, 5);
	const std::string s3 = "CTCCTGACCTCAGGTGATCTGCTCATCTTGGCCTCTCAAAGTGCTGGGATGACAGGTGTGAGCCACCGTGCCTGGCTTTTTTTTTTTTTTCTTTCCTATTTCCCCATTGGGTCCGATTCCCCTGCCTCAGCCTCCCGAGCAGCTGGGACTACAGGCACCCGCCACCACGCCCGGCTAATTTTTGTATTTTTAGTAGAGACGGGGTTTCACCATGTTGGCCAGGCTGGTCGCGAACTCCTGACCTCAGGTGATCCGCCCGCCTCGGCCTCCCAAAGAGTGTGTTCATGTCTCACCAGCACAAATATGTAATCAACCGTGTCCCCTTTGCAGAGTAAGGTGGCCATGAGGTTTCAGGTATATTATCGGATCCTTCTATTCCATTTTAATGACCTCCCATTATTTCTGCAGAGAAGGTGAAAGTTCCCCTTCAACATCCTTCATCACAGCCATCGCCCCCCAACTCCACTCGCCACTGCCACCCAGACGTCAAATCCTAGAGAGACTTCAGATGGCCGGGGGAGGCCTTCCTGGCTTTCTTTGTCACTCCACGCCGAATCCGAGTGTCACACGTGGCCCACCTGTGAAGTTCCTGGCCGGTGCTGATATAAATGGTAACAGGGCCGATGGCCTGCGATTGTCAGAGATTGATTTAGACAGGTGACTTCTTCTGCTTTGTAAGGTAATGAGACCGAACTTTCAGACTAACTTTCCAAACTCTGGGGGTGACGGCTTCAGTGCAGTCAGTTGACGTAAAAGTTTATCGGCAGGCAGCGGTGGCTCACGCCTGGAATCCCACAACTTTAGTAGAGACAGGGTTTCACCGTGTTGGCCAGGCTGGTCTCAAACGCCTGACCTAACGTGATCCACCTGTCTTGGCCTCCCAGAGTGCTGGAATTACGACCGTGAGCCACCACGCCTGGCCCAACAAAACTTATTTAATTTTTATTTTATTTGTTTATTTTTTGAGACAGAGTCCCACTCTGTCGCCTAGGCTGGAGTACAGTGGCGCAATCTCGGCTCACT";
	// ��NC_000024.10 Homo sapiens chromosome Y, GRCh38.p14�Ͻ�ȡһ��
	OrfScanning(std::cout, s3, 30);
	return 0;
}


struct Orf
{
	int start, end;
	Orf(int s, int e):start(s), end(e){}
	int len() { return abs(start - end) + 1;}
};

// ����Ϊ���������ɨ��Ĵ���ʶ��ΪORF�ĳ�����ֵ
void OrfScanning(std::ostream& os, const std::string &seq, int len)
{
	std::list<Orf> OrfLst; // ����ɨ�赽��ORF������
	std::string coden, re_coden; // ˫���ȡ������������

	os << "ORF scanning (len >= " << len << ")\nseq: " << seq << '\n';

	for (int i = 0; i < 3; i++) // 3��ORF
	{
		for (int end = i; end + 2 < seq.size(); end += 3) // ɨ��һ��ORF������������
		{
			coden = seq.substr(end, 3);
			re_coden = std::string(coden.rbegin(), coden.rend()); // �����ȡ

			if (IsTermination(coden)) // ����ҵ���ֹ������
			{
				for (int start = end - 3; start >= 0; start -= 3) // �Ӻ���ǰѰ����ʼ������
				{
					coden = seq.substr(start, 3);
					if (IsTermination(coden)) // ����ҵ���ֹ������
						break;
					if (IsInitiation(coden)) // ����ҵ���ʼ������
					{
						OrfLst.push_front(Orf(start + 1, end + 3)); // ���ҵ���ORF��λ�ñ��浽�����У��±��1��ʼ
						break;
					}   
				}
			}

			if (IsTermination(re_coden)) // ����ҵ��������ֹ������
			{
				for (int start = end + 3; start + 2 < seq.size(); start += 3) // �Ӻ���ǰѰ����ʼ������
				{
					coden = seq.substr(start, 3);
					re_coden = std::string(coden.rbegin(), coden.rend()); // �����ȡ
					if (IsTermination(re_coden)) // ����ҵ���ֹ������
						break;
					if (IsInitiation(re_coden)) // ����ҵ���ʼ������
					{
						OrfLst.push_front(Orf(start + 3, end + 1)); // ���ҵ���ORF��λ�ñ��浽�����У��±��1��ʼ
						break;
					}
				}
			}
		} // ɨ��һ��ORF������������
	} // 3��ORF

	while (!OrfLst.empty()) // ����������ȴ�����ֵ��ORF����
	{
		if (OrfLst.front().len() >= len)
		{
			os << "start:" << OrfLst.front().start << "\tend:" << OrfLst.front().end << "\tlen:" << OrfLst.front().len() << '\n';
			os << seq.substr((OrfLst.front().start < OrfLst.front().end ? OrfLst.front().start : OrfLst.front().end) - 1, OrfLst.front().len()) << '\n';
		}
		OrfLst.pop_front();
	}
	os << "ORF scanning is complete!\n" << std::endl;
}

bool IsInitiation(std::string coden)
{
	assert(coden.size() == 3);
	std::transform(coden.begin(), coden.end(), coden.begin(), ::toupper); //ת��д
	if (coden.compare(std::string("ATG")) == 0) //��ʼ�����Ӷ�ӦDNA��ATG
		return true;
	else
		return false;
}

bool IsTermination(std::string coden)
{
	assert(coden.size() == 3);
	std::transform(coden.begin(), coden.end(), coden.begin(), ::toupper); //ת��д
	// ��ֹ�����Ӷ�ӦDNA��TAG,TAA,TGA
	if ((coden.compare(std::string("TAG")) && coden.compare(std::string("TAA")) && coden.compare(std::string("TGA"))) == false)
		return true;
	else
		return false;
}
