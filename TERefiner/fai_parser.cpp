#include<fstream>
#include"fai_parser.h"
#include<iostream>

using namespace std;

FaiParser::FaiParser()
{
	this->path="";
}
	
FaiParser::FaiParser(std::string path)
{
	this->path=path;
}

void FaiParser::setPath(std::string path)
{
	this->path=path;
}

/*
Description
	Parse .fasta.fai file, and returen chromosome start position and number of chars each line.
*/
void FaiParser::parseFai()
{
	ifstream fin;
	fin.open(path.c_str(), ifstream::in);

	string chrom_name;
	int length;
	unsigned long start_pos;
	int size_symbol;
	int size_ascii;
	int cnt=0;
	vchroms.clear();
	while(fin>>chrom_name>>length>>start_pos>>size_symbol>>size_ascii)
	{

		//cout<<chrom_name<<" "<<length<<" "<<start_pos<<" "<<size_symbol<<" "<<size_ascii<<endl;//////////////////////////////////////////////////////

		ChromInfo ci;
		ci.id=cnt;
		ci.cname=chrom_name;
		ci.length=length;
		ci.startpos=start_pos;
		ci.size_chars=size_symbol;
		ci.size_ascii=size_ascii;
		vchroms.push_back(ci);
		
		cnt++;
	}
	fin.close();
}

/*
Description:
	Parse .fasta.fai file, and returen chromosome length.
Input:
	chrom name
Output:
	length of chrom
*/
int FaiParser::getChromLen(std::string chrom)
{	
	int size=vchroms.size();
	for(int i=0;i<size;i++)
	{
		if(vchroms[i].cname==chrom)
		{
			return vchroms[i].length;
		}
	}
	return 0;
}

string FaiParser::getChromName(int chrom_id)
{
	return this->vchroms[chrom_id].cname;
}