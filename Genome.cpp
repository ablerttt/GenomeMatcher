#include "provided.h"
#include <string>
#include <vector>
#include <iostream>
#include <istream>
using namespace std;

class GenomeImpl
{
public:
	GenomeImpl(const string& nm, const string& sequence);
	static bool load(istream& genomeSource, vector<Genome>& genomes);
	int length() const;
	string name() const;
	bool extract(int position, int length, string& fragment) const;
private:
	string m_name;
	string m_sequence;
	int m_length;
};

GenomeImpl::GenomeImpl(const string& nm, const string& sequence)
{
	m_name = nm;
	m_sequence = sequence;
	m_length = m_sequence.size();
}

bool GenomeImpl::load(istream& genomeSource, vector<Genome>& genomes)
{
	genomes.clear();
	if (!genomeSource) //invalid istream 
		return false;

	string currentLine;
	string currentGenome;
	string currentName;

	while (getline(genomeSource, currentLine))
	{
		if (currentLine[0] == '>') //indicates that we're looking at the name
		{
			if (currentGenome.size() > 0 && currentName.size() > 0)
			{
				Genome g(currentName, currentGenome);
				genomes.push_back(g);

				currentName = currentLine.substr(1);
				if (currentName.size() == 0)
					return false;
				currentGenome = "";
			}
			else if (currentGenome.size() > 0 && currentName.size() == 0) //genome with no name
				return false;
			else if (currentGenome.size() == 0 && currentName.size() > 0) //name after a name
				return false;
			else //current genome and current name are empty
			{
				currentName = currentLine.substr(1);
				//name line with a greater than character but no other characters
				if (currentName.size() == 0)
					return false;
			}
		}
		else
		{
			//meaning that this is another line for the genome
			for (char ch : currentLine)
			{
				switch (ch)
				{
					case 'a': case 'A': currentGenome.append("A"); break;
					case 't': case 'T': currentGenome.append("T"); break;
					case 'c': case 'C': currentGenome.append("C"); break;
					case 'g': case 'G': currentGenome.append("G"); break;
					case 'n': case 'N': currentGenome.append("N"); break;
					default: 			return false; 		 	   break;
				}
			}
		}
	}
	if (currentName != "" && currentGenome != "")
	{				
		Genome g(currentName, currentGenome);
		genomes.push_back(g);
	}
	return true;
}

int GenomeImpl::length() const
{
	return m_length;
}

string GenomeImpl::name() const
{
	return m_name;
}

bool GenomeImpl::extract(int position, int length, string& fragment) const
{
	if (position + length > m_length)
		return false;
	fragment = "";

	fragment = m_sequence.substr(position, length);
	return true;
}

//******************** Genome functions ************************************

// These functions simply delegate to GenomeImpl's functions.
// You probably don't want to change any of this code.

Genome::Genome(const string& nm, const string& sequence)
{
	m_impl = new GenomeImpl(nm, sequence);
}

Genome::~Genome()
{
	delete m_impl;
}

Genome::Genome(const Genome& other)
{
	m_impl = new GenomeImpl(*other.m_impl);
}

Genome& Genome::operator=(const Genome& rhs)
{
	GenomeImpl* newImpl = new GenomeImpl(*rhs.m_impl);
	delete m_impl;
	m_impl = newImpl;
	return *this;
}

bool Genome::load(istream& genomeSource, vector<Genome>& genomes)
{
	return GenomeImpl::load(genomeSource, genomes);
}

int Genome::length() const
{
	return m_impl->length();
}

string Genome::name() const
{
	return m_impl->name();
}

bool Genome::extract(int position, int length, string& fragment) const
{
	return m_impl->extract(position, length, fragment);
}
