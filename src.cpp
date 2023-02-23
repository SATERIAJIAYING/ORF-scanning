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
	// 从NC_000024.10 Homo sapiens chromosome Y, GRCh38.p14上截取一段
	OrfScanning(std::cout, s3, 30);
	return 0;
}


struct Orf
{
	int start, end;
	Orf(int s, int e):start(s), end(e){}
	int len() { return abs(start - end) + 1;}
};

// 参数为输出流、带扫描的串、识别为ORF的长度阈值
void OrfScanning(std::ostream& os, const std::string &seq, int len)
{
	std::list<Orf> OrfLst; // 保存扫描到的ORF的链表
	std::string coden, re_coden; // 双向读取的三联密码子

	os << "ORF scanning (len >= " << len << ")\nseq: " << seq << '\n';

	for (int i = 0; i < 3; i++) // 3个ORF
	{
		for (int end = i; end + 2 < seq.size(); end += 3) // 扫描一个ORF的三联密码子
		{
			coden = seq.substr(end, 3);
			re_coden = std::string(coden.rbegin(), coden.rend()); // 反向读取

			if (IsTermination(coden)) // 如果找到终止密码子
			{
				for (int start = end - 3; start >= 0; start -= 3) // 从后往前寻找起始密码子
				{
					coden = seq.substr(start, 3);
					if (IsTermination(coden)) // 如果找到终止密码子
						break;
					if (IsInitiation(coden)) // 如果找到起始密码子
					{
						OrfLst.push_front(Orf(start + 1, end + 3)); // 将找到的ORF的位置保存到链表中，下标从1开始
						break;
					}   
				}
			}

			if (IsTermination(re_coden)) // 如果找到反向的终止密码子
			{
				for (int start = end + 3; start + 2 < seq.size(); start += 3) // 从后往前寻找起始密码子
				{
					coden = seq.substr(start, 3);
					re_coden = std::string(coden.rbegin(), coden.rend()); // 反向读取
					if (IsTermination(re_coden)) // 如果找到终止密码子
						break;
					if (IsInitiation(re_coden)) // 如果找到起始密码子
					{
						OrfLst.push_front(Orf(start + 3, end + 1)); // 将找到的ORF的位置保存到链表中，下标从1开始
						break;
					}
				}
			}
		} // 扫描一个ORF的三联密码子
	} // 3个ORF

	while (!OrfLst.empty()) // 依次输出长度大于阈值的ORF到流
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
	std::transform(coden.begin(), coden.end(), coden.begin(), ::toupper); //转大写
	if (coden.compare(std::string("ATG")) == 0) //起始密码子对应DNA：ATG
		return true;
	else
		return false;
}

bool IsTermination(std::string coden)
{
	assert(coden.size() == 3);
	std::transform(coden.begin(), coden.end(), coden.begin(), ::toupper); //转大写
	// 终止密码子对应DNA：TAG,TAA,TGA
	if ((coden.compare(std::string("TAG")) && coden.compare(std::string("TAA")) && coden.compare(std::string("TGA"))) == false)
		return true;
	else
		return false;
}
