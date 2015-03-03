#include"fasta_parser.h"
#include<fstream>
#include"public_parameters.h"


using namespace std;

//const int MAX_LEN_EACH_LINE_FAI=300;

FastaParser::FastaParser()
{
	this->vid_name.clear();
}

FastaParser::FastaParser(std::string path)
{
	this->path=path;
	this->vid_name.clear();
}

void FastaParser::setPath(std::string path)
{
	this->path=path; 
}

/*
Description: return the specified sub-sequence.
*/
string FastaParser::parseFasta(std::string chrom,int pos, int length)
{	
	int chrom_start=0;//start position of chromosome, including all ascii chars, like '\0''\n'.
	int size_each_line=0;//number of symbols, like without  '\0\n'.
	int size_ascii=0;//number of all ascii chars, like with  '\0\n'.
	
	parseFai(chrom, chrom_start, size_each_line, size_ascii);//get chrom start position from fasta.fai file.
	if(chrom_start==0)
	{
		return "";
	}

	ifstream fin;
	fin.open(path.c_str(), ifstream::in);//open fasta file


	int start_part_len=size_each_line-size_ascii;
	int lines=0;
	int left=pos%size_each_line;
	if(pos>size_each_line)
	{
		fin.seekg(chrom_start);//move to chrom_start position
		string start_part;
		fin>>start_part;
		start_part_len=start_part.length();

		lines=(pos-start_part_len)/size_each_line;
		left=(pos-start_part_len)%size_each_line;

	}
	int newpos=(start_part_len+(size_ascii-size_each_line)) + lines*size_ascii + left;
	
	fin.seekg(newpos+chrom_start);

	string sub_ref="";
	int cnt=0;
	while(true)
	{
		string temp;
		fin>>temp;
		cnt+=temp.length();
		if(cnt>=length)
		{
			sub_ref+=temp.substr(0,(size_each_line-(cnt-length)));
			break;
		}
		else
		{
			sub_ref+=temp;
		}
	}
	fin.close();

	return sub_ref;
}

//void FastaParser::mapChromIDName()//map chrom ID with Name, by using XX.fasta.fai file.
//{
//	ifstream fin;
//	std::string fai_path=path+".fai";
//	fin.open(fai_path.c_str(), ifstream::in);
//
//	int id=0;
//	ofstream fout;
//	fout.open(CHROM_ID_NAME.c_str(),ofstream::out);
//	std::string chrom_name;
//	while(fin>>chrom_name)
//	{
//		char buffer[MAX_LEN_EACH_LINE_FAI];
//		fin.getline(buffer,MAX_LEN_EACH_LINE_FAI);
//		
//		fout<<id<<" "<<chrom_name<<endl;
//		id++;
//	}
//	fout.close();
//	fin.close();
//}
//
//void FastaParser::loadChromIDName()
//{
//	ifstream fin;
//	fin.open(CHROM_ID_NAME.c_str(), ifstream::in);
//	int id;
//	string name;
//	vid_name.clear();
//	while(fin>>id>>name)
//	{
//		this->vid_name.push_back(std::make_pair(id,name));
//	}
//	fin.close();
//}

//------------private functions------------------------------------------------------
/*
Description:
	Parse .fasta.fai file, and returen chromosome start position and number of chars each line.
*/
void FastaParser::parseFai(std::string chrom, int& start, int& size_aline, int& size_all_aline)
{
	ifstream fin;
	std::string fai_path=path+".fai";
	fin.open(fai_path.c_str(), ifstream::in);
	
	string chrom_name;
	int length;
	int start_pos;
	int size_symbol;
	int size_ascii;
	
	while(fin>>chrom_name>>length>>start_pos>>size_symbol>>size_ascii)
	{
		if(chrom_name==chrom)
		{
			start=start_pos;
			size_aline=size_symbol;
			size_all_aline=size_ascii;
			break;
		}
	}
	fin.close();
}