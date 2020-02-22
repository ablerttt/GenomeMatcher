#ifndef TRIE_INCLUDED
#define TRIE_INCLUDED

#include <string>
#include <algorithm>
#include <vector>

template<typename ValueType>
class Trie
{
public:
	Trie();
	~Trie();
	void reset();
	void insert(const std::string& key, const ValueType& value);
	std::vector<ValueType> find(const std::string& key, bool exactMatchOnly) const;

	// C++11 syntax for preventing copying and assignment
	Trie(const Trie&) = delete;
	Trie& operator=(const Trie&) = delete;
private:
	void recursivelyReset();
	std::vector<ValueType> findWithError(const std::string& key, bool oneError, bool first) const;
	int vContains(std::vector<char> input, char c) const;
	struct Node
	{
		std::vector<std::vector<ValueType>> values;
		std::vector<Trie*> children;
		std::vector<char> childrenValue;
	};
	Node* root;
};

#endif // TRIE_INCLUDED

template<typename ValueType>
Trie<ValueType>::Trie()
{
	root = new Node;
}

template<typename ValueType>
Trie<ValueType>::~Trie()
{
	for (auto p : root->children)
		delete p;
	delete root;
}

template<typename ValueType>
void Trie<ValueType>::reset()
{
	recursivelyReset();
	root = new Node();
}

template<typename ValueType>
void Trie<ValueType>::recursivelyReset()
{
	if (root == nullptr)
	{
		return;
	}

	for (int i = 0; i < static_cast<int>(root->children.size()); i++)
	{
		root->children[i]->recursivelyReset();
	}

	delete root;
}

template<typename ValueType>
void Trie<ValueType>::insert(const std::string & key, const ValueType & value)
{
	int contains = vContains(root->childrenValue, key[0]);

	if (contains == -1)
	{
		root->childrenValue.push_back(key[0]);
		root->values.push_back(std::vector<ValueType>());
		root->children.push_back(new Trie());
	}

	if (key.size() == 1)
	{
		if (contains > -1)
			root->values[contains].push_back(value);
		else
			root->values[root->values.size() - 1].push_back(value);
		return;
	}

	if (contains > -1)
		root->children[contains]->insert(key.substr(1), value);
	else
		root->children[root->children.size() - 1]->insert(key.substr(1), value);
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::find(const std::string & key, bool exactMatchOnly) const
{
	if (!exactMatchOnly)
		return findWithError(key, true, true);

	int path = vContains(root->childrenValue, key[0]);
	if (path == -1)
		return std::vector<ValueType>();

	if (key.size() == 1)
		return root->values[path];

	return root->children[path]->find(key.substr(1), true);
}

template<typename ValueType>
std::vector<ValueType> Trie<ValueType>::findWithError(const std::string & key, bool oneError, bool first) const
{
	int path = vContains(root->childrenValue, key[0]);
	if (first && path == -1)
		return std::vector<ValueType>();
	if (!oneError && path == -1)
		return std::vector<ValueType>();
	if (key.size() == 1 && path == -1)
		return std::vector<ValueType>();

	if (key.size() == 1 && oneError)
	{
		std::vector<ValueType> lastOnes;
		for (int i = 0; i < static_cast<int>(root->values.size()); i++)
		{
			for (int j = 0; j < static_cast<int>(root->values[i].size()); j++)
				lastOnes.push_back(root->values[i][j]);
		}

		return lastOnes;
	}
	else if (key.size() == 1 && !oneError)
		return root->values[path];

	std::vector<ValueType> best;
	if (first)
		best = root->children[path]->findWithError(key.substr(1), true, false);
	else if (!oneError)
		best = root->children[path]->findWithError(key.substr(1), false, false);
	else if (oneError)
	{
		for (int i = 0; i < static_cast<int>(root->children.size()); i++)
		{
			std::vector<ValueType> temp;
			if (i == path)
				temp = root->children[i]->findWithError(key.substr(1), true, false);
			else
				temp = root->children[i]->findWithError(key.substr(1), false, false);

			for (int j = 0; j < static_cast<int>(temp.size()); j++)
				best.push_back(temp[j]);
		}
	}
	return best;
}

template<typename ValueType>
int Trie<ValueType>::vContains(std::vector<char> input, char c) const
{
	int count = 0;
	for (std::vector<char>::iterator it = input.begin(); it != input.end(); it++)
	{
		if ((*it) == c)
			return count;

		count++;
	}
	return -1;
}

