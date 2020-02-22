#include "provided.h"
#include "Trie.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <utility>
#include <algorithm>
using namespace std;

class GenomeMatcherImpl
{
public:
	GenomeMatcherImpl(int minSearchLength);
	void addGenome(const Genome& genome);
	int minimumSearchLength() const;
	bool findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const;
	bool findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const;
private:
	int m_minLength;
	Trie<pair<int, int>> positions;
	vector<Genome*>* genomes;
	int genomeLoc(string name) const;
	int genomePosition = 0;
};

bool compare(const GenomeMatch& a, const GenomeMatch& b);

GenomeMatcherImpl::GenomeMatcherImpl(int minSearchLength)
{
	m_minLength = minSearchLength;
	genomes = new vector<Genome*>();
}

void GenomeMatcherImpl::addGenome(const Genome& genome)
{
	genomes->push_back(new Genome(genome));
	string temp;
	for (int i = 0; i < genome.length(); i++)
	{
		if (genome.extract(i, m_minLength, temp)) {
			pair<int, int> insertedPair(genomePosition, i);
			positions.insert(temp, insertedPair);
		}
	}
	genomePosition++;
}

int GenomeMatcherImpl::minimumSearchLength() const
{
	return m_minLength;
}

bool GenomeMatcherImpl::findGenomesWithThisDNA(const string& fragment, 
	int minimumLength, 
	bool exactMatchOnly, 
	vector<DNAMatch>& matches) const
{
	if (fragment.size() < minimumLength) 
		return false;

	if (minimumLength < m_minLength)
		return false;

	vector<pair<int, int>> duplicates = positions.find(fragment.substr(0, m_minLength), exactMatchOnly); //first is genomeloc, second is position
	vector<pair<int, int>> maxSizes(genomes->size()); //first is position, second is length

	
	for (vector<pair<int, int>>::iterator it = maxSizes.begin(); it != maxSizes.end(); it++)
	{
		it->first = - 1;
		it->second = minimumLength;
	}

	for (vector<pair<int, int>>::iterator it = duplicates.begin(); it != duplicates.end(); it++)
	{
		int genomeLoc = it->first;
		string embeddedFragment;
		for (int longestPossible = fragment.length(); longestPossible >= minimumLength; longestPossible--)
		{
			if ((*genomes)[genomeLoc]->extract(it->second, longestPossible, embeddedFragment))
				break;
		}

		int length = embeddedFragment.size();
		int maxErrors = 0;
		if (!exactMatchOnly) maxErrors++;

		for (int i = 0; i < length; i++)
		{
			if (i == 0 && fragment[0] != embeddedFragment[0])
				continue;
			if (exactMatchOnly && fragment[i] != embeddedFragment[i])
			{
				length = i;
				break;
			}
			else if (!exactMatchOnly && fragment[i] != embeddedFragment[i])
			{
				maxErrors--;
			}
			if (maxErrors < 0)
			{
				length = i;
				break;
			}
		}
		if (length >= maxSizes[genomeLoc].second)
		{
			if (maxSizes[genomeLoc].first >= 0 && it->second > maxSizes[genomeLoc].first)
				continue;
			maxSizes[genomeLoc].first = it->second;
			maxSizes[genomeLoc].second = length;
		}
	}	 

	vector<bool> taken(genomes->size());
	for (vector<bool>::iterator it = taken.begin(); it != taken.end(); it++)
		*it = false;

	for (vector<pair<int, int>>::iterator it = duplicates.begin(); it != duplicates.end(); it++)
	{
		int genomeLoc = it->first;
		if (it->second == maxSizes[genomeLoc].first && taken[genomeLoc] == false)
		{
			DNAMatch d;
			d.genomeName = (*genomes)[genomeLoc]->name();
			d.length = maxSizes[genomeLoc].second;
			d.position = maxSizes[genomeLoc].first;
			matches.push_back(d);
			taken[genomeLoc] = true;
		}
	}

	
	return matches.size() > 0;
}

bool GenomeMatcherImpl::findRelatedGenomes(const Genome& query, 
	int fragmentMatchLength, 
	bool exactMatchOnly, 
	double matchPercentThreshold, 
	vector<GenomeMatch>& results) const
{
	if (fragmentMatchLength < m_minLength)
		return false;

	vector<int> matches(genomes->size());
	for (vector<int>::iterator it = matches.begin(); it != matches.end(); it++) { *it = 0; }

	int numFragments = query.length() / fragmentMatchLength;

	for (int pos = 0; pos < query.length(); pos += fragmentMatchLength)
	{
		string fragment;
		if (!query.extract(pos, fragmentMatchLength, fragment)) continue;

		vector<DNAMatch> match;
		if (!findGenomesWithThisDNA(fragment, fragmentMatchLength, exactMatchOnly, match)) continue;

		for (vector<DNAMatch>::iterator it = match.begin(); it != match.end(); it++)
		{
			int loc = genomeLoc(it->genomeName);
			matches[loc]++;
		}
	}

	for (int i = 0; i < genomes->size(); i++)
	{
		if (matches[i] == 0) continue;
		double percentage = 100.0 * matches[i] / numFragments;
		if (percentage >= matchPercentThreshold)
		{
			GenomeMatch gm;
			gm.genomeName = (*genomes)[i]->name();
			gm.percentMatch = percentage;
			results.push_back(gm);
		}
	}

	sort(results.begin(), results.end(), compare);
	return results.size() > 0;
}


int GenomeMatcherImpl::genomeLoc(string name) const
{
	for (int i = 0; i < genomes->size(); i++)
		if ((*genomes)[i]->name() == name) 
			return i;
	return -1;
}

bool compare(const GenomeMatch& a, const GenomeMatch& b)
{
	if (a.percentMatch > b.percentMatch)
		return true;
	else if (a.percentMatch < b.percentMatch)
		return false;
	else
		return (a.genomeName < b.genomeName);
}

//******************** GenomeMatcher functions ********************************

// These functions simply delegate to GenomeMatcherImpl's functions.
// You probably don't want to change any of this code.

GenomeMatcher::GenomeMatcher(int minSearchLength)
{
	m_impl = new GenomeMatcherImpl(minSearchLength);
}

GenomeMatcher::~GenomeMatcher()
{
	delete m_impl;
}

void GenomeMatcher::addGenome(const Genome& genome)
{
	m_impl->addGenome(genome);
}

int GenomeMatcher::minimumSearchLength() const
{
	return m_impl->minimumSearchLength();
}

bool GenomeMatcher::findGenomesWithThisDNA(const string& fragment, int minimumLength, bool exactMatchOnly, vector<DNAMatch>& matches) const
{
	return m_impl->findGenomesWithThisDNA(fragment, minimumLength, exactMatchOnly, matches);
}

bool GenomeMatcher::findRelatedGenomes(const Genome& query, int fragmentMatchLength, bool exactMatchOnly, double matchPercentThreshold, vector<GenomeMatch>& results) const
{
	return m_impl->findRelatedGenomes(query, fragmentMatchLength, exactMatchOnly, matchPercentThreshold, results);
}
