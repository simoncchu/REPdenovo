#ifndef _H_RESRAP_ATSAF_
#define _H_RESRAP_ATSAF_

#include<string>
#include<vector>
#include<utility>
#include"fai_parser.h"

/*
Description:
	Parse given fasta file by analyzing "xx.fai" file.
*/
class FastaParser
{
public:
	FastaParser();
	FastaParser(std::string path);

public:
	std::string parseFasta(std::string chrom,int pos, int length);
	//void mapChromIDName();//map chrom ID with Name, by using XX.fasta.fai file 
	//void loadChromIDName();//load chrom name and id into memory

public:
	void setPath(std::string path);

private:
	void parseFai(std::string chrom, int& start, int& size_each_line, int& size_ascii);

public:
	std::vector<std::pair<int,std::string> > vid_name;//save relationship between chromosome name and id.

private:
	std::string path;
};

#endif
