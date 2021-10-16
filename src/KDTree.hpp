#ifndef SRC_KDTREE_HPP_
#define SRC_KDTREE_HPP_

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>
#include "Point.hpp"

template <typename value_type>
struct KDTreeNode
{
    value_type nodeValue;
    KDTreeNode<value_type>* arrNod[2];
    KDTreeNode(const value_type& value)
    {
        nodeValue = value;
        arrNod[0] = arrNod[1] = 0;
    }
    KDTreeNode(value_type value, KDTreeNode<value_type>* left, KDTreeNode<value_type>* right)
    {
        nodeValue = value;
        arrNod[0] = left;
        arrNod[1] = right;
    }
};

template <size_t N, typename ElemType>
class KDTree {
public:
    typedef std::pair<Point<N>, ElemType> value_type;

    KDTree();
    ~KDTree();

    KDTree(const KDTree& rhs);
    KDTree& operator=(const KDTree& rhs);

    size_t dimension() const;

    size_t size() const;
    bool empty() const;

    bool find(const Point<N>& pt, KDTreeNode<value_type>**& p) const;

    bool contains(const Point<N>& pt) const;

    void insert(const Point<N>& pt, const ElemType& value);

    ElemType& operator[](const Point<N>& pt);

    ElemType& at(const Point<N>& pt);
    const ElemType& at(const Point<N>& pt) const;

    void knn_aux(Point<N> key, KDTreeNode<value_type>* currentNode, int nivel, std::vector<std::pair<ElemType, double> >& cTipo) const;
    ElemType knn_value(const Point<N>& key, size_t k) const;

    std::vector<ElemType> knn_query(const Point<N>& key, size_t k) const;

private:
    size_t dimension_;
    size_t size_;
    KDTreeNode<value_type>* root = nullptr;
};

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree()
{
    dimension_ = N;
    size_ = 0;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::~KDTree()
{
    eraseKDTree(root);
}

template <typename value_type>
void eraseKDTree(KDTreeNode<value_type>* node)
{
    if (node != nullptr) {
        eraseKDTree((node->arrNod)[1]);
        eraseKDTree((node->arrNod)[0]);
        delete node;
    }
}

template <typename value_type>
KDTreeNode<value_type>* copyNodes(const KDTreeNode<value_type>* node)
{
    if (node != nullptr)
    {
        KDTreeNode<value_type>* nodeCopy = new KDTreeNode<value_type>(node->nValue, copyNodes((node->arrNod)[0]), copyNodes((node->arrNod)[1]));
        return nodeCopy;
    }
    return nullptr;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>::KDTree(const KDTree& rhs)
{
    root = copyNodes(rhs.root);
    dimension_ = rhs.dimension_;
    size_ = rhs.size_;
}

template <size_t N, typename ElemType>
KDTree<N, ElemType>& KDTree<N, ElemType>::operator=(const KDTree& rhs)
{
    root = copyNodes(rhs.root);
    dimension_ = rhs.dimension_;
    size_ = rhs.size_;
    return *this;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::dimension() const
{
    return N;
}

template <size_t N, typename ElemType>
size_t KDTree<N, ElemType>::size() const
{
    return size_;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::empty() const
{
    return !root;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::find(const Point<N>& pt, KDTreeNode<value_type>**& p) const
{
    size_t i = 0;
    for (p = const_cast<KDTreeNode<value_type>**>(&root);*p &&((*p)->nValue).first!=pt;p=&((*p)->arrNod[pt[i % N]>(((*p)->nValue).first)[i%N]])){i++;}
    return *p != 0;
}

template <size_t N, typename ElemType>
bool KDTree<N, ElemType>::contains(const Point<N>& pt) const
{
    KDTreeNode<value_type>** p;
    if (!find(pt, p))
    {
        return false;
    }
    return true;
}

template <size_t N, typename ElemType>
void KDTree<N, ElemType>::insert(const Point<N>& pt, const ElemType& value)
{
    KDTreeNode<value_type>** p;
    if (!find(pt, p))
    {
        value_type aux;
        aux.first = pt;
        aux.second = value;
        *p = new KDTreeNode<value_type>(aux);
        size_ += 1;
    }
    ((*p)->nValue).second = value;
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::operator[](const Point<N>& pt)
{
    KDTreeNode<value_type>** p;
    if (!find(pt, p))
    {
        value_type aux;
        aux.first = pt;
        aux.second = size_;
        *p = new KDTreeNode<value_type>(aux);
        size_ += 1;
        return ((*p)->nValue).second;
    }
    return ((*p)->nValue).second;
}

template <size_t N, typename ElemType>
ElemType& KDTree<N, ElemType>::at(const Point<N>& pt)
{
    KDTreeNode<value_type>** p;
    if (find(pt, p))
    {
        return ((*p)->nValue).second;
    }
    throw std::out_of_range("");
}

template <size_t N, typename ElemType>
const ElemType& KDTree<N, ElemType>::at(const Point<N>& Pt) const
{
    KDTreeNode<value_type>** P;
    if (find(Pt, P))
    {
        return ((*P)->nValue).second;
    }
    throw std::out_of_range("");
}


template <size_t N, typename ElemType>
void KDTree<N, ElemType>::knn_aux(Point<N> key, KDTreeNode<value_type>* cNode, int nivel, std::vector<std::pair<ElemType, double> >& cTipo) const
{
    if (cNode == nullptr)
    {
        return;
    }
    double d = distance((cNode->nValue).first, key);
    cTipo.push_back(std::make_pair((cNode->nValue).second, d));

    knn_aux(key, (cNode->arrNod)[0], ++nivel, cTipo);
    knn_aux(key, (cNode->arrNod)[1], ++nivel, cTipo);
}
template <size_t N, typename ElemType>
ElemType KDTree<N, ElemType>::knn_value(const Point<N>& key, size_t k) const
{
    if (k > size_)
    {
        k = size_;
    }
    std::vector<ElemType> values(k);
    values = knn_query(key, k);
    std::vector<ElemType> ElemV;
    ElemV.push_back(values[0]);
    for (size_t i = 1; i < k; i++)
    {
        bool bol = false;
        for (size_t j = 0; j < ElemV.size(); j++)
        {
            if (values[i] == ElemV[j])
            {
                bol = true;
            }
        }
        if (!ex)
        {
            ElemV.push_back(values[i]);
        }
    }
    std::vector<std::pair<ElemType, size_t> > cont;
    for (size_t i = 0; i < ElemV.size(); i++)
    {
        cont.push_back(std::make_pair(ElemV[i], 0));
    }
    for (size_t i = 0; i < k; i++)
    {
        for (size_t j = 0; j < cont.size(); j++)
        {
            if (values[i] == cont[j].first)
            {
                cont[j].second += 1;
            }
        }
    }
    std::pair<ElemType, size_t> maxi = cont[0];
    for (size_t i = 1; i < cont.size(); i++)
    {
        if (maxi.second < cont[i].second)
        {
            maxi = cont[i];
        }
    }
    return maxi.first;
}

template <size_t N, typename ElemType>
std::vector<ElemType> KDTree<N, ElemType>::knn_query(const Point<N>& key, size_t k) const
{
    std::vector<std::pair<ElemType, double> > cTipo;
    KDTreeNode<value_type>* cNode = root;
    size_t nivel = 0;

    knn_aux(key, cNode, nivel, cTipo);
    std::sort(cTipo.begin(), cTipo.end(), [](const auto& x, const auto& y) { return x.second < y.second; });

    std::vector<ElemType> val;
    for (size_t i = 0; i < k; i++) {
        val.push_back(cTipo[i].first);
    }
    return val;
}

#endif  // SRC_KDTREE_HPP_
