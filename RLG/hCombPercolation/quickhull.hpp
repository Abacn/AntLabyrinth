/* Quickhull algorithm implementation
 *
 * Copyright (c) 2014-2015, Anatoliy V. Tomilov
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following condition is met:
 * Redistributions of source code must retain the above copyright notice, this condition and the following disclaimer.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * Modefied by Yi Hu (2020), compatible with c++11, fix dimension given by template parameter,
 * reduce memeory usage by replacing std::vector to std::array, std::list to std::forward_list, and etc
 */
#pragma once

#include <type_traits>
#include <array>
#include <stack>
#include <vector>
#include <deque>
#include <list>
#include <forward_list>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <iterator>
#include <memory>
#include <algorithm>
#include <numeric>
#include <utility>
#include <functional>

#include <cstdint>
#include <cmath>
#include <cassert>

template< int D, typename point_iterator,
typename value_type = typename std::decay< decltype(*std::begin(std::declval< typename std::iterator_traits< point_iterator >::value_type >())) >::type >
struct quick_hull
{
  
  static_assert(std::is_base_of< std::forward_iterator_tag, typename std::iterator_traits< point_iterator >::iterator_category >::value,
                "multipass guarantee required");
  
  using small_int_type = std::int32_t;
  using small_uint_type = std::uint32_t;
  using size_type = std::size_t;
  //size_type const dimension_;
  value_type const & eps;
  
  value_type const zero = value_type(0);
  value_type const one = value_type(1);
  
  using vector = std::array< value_type, D >;
  
  private :
  
  using vrow = value_type *;
  using crow = value_type const *;
  using matrix = std::array< vrow, D >;
  
  value_type storage_[ D * D * 2 + D];
  vrow inner_point_;
  matrix matrix_;
  matrix det_matrix_;
  matrix shadow_matrix_;
  
  public :
  
  quick_hull(value_type const &&) = delete; // bind eps to lvalue only
  
  quick_hull(value_type const & _eps)
  : eps(_eps)
  , inner_point_(storage_)
  {
    assert(1 < D);
    assert(!(eps < zero));
    for (vrow & row_ : matrix_) {
      row_ = inner_point_;
      inner_point_ += D;
    }
    for (vrow & row_ : shadow_matrix_) {
      row_ = inner_point_;
      inner_point_ += D;
    }
    //assert(inner_point_ + D == &storage_.back() + 1);
  }
  
  using point_array = std::array< point_iterator, D >;
  using point_list  = std::list< point_iterator >;
  using point_forward_list = std::forward_list< point_iterator >;
  using point_deque = std::deque< point_iterator >;
  
  using facet_array = std::array< small_uint_type, D >;
  using facet_container = std::vector< small_uint_type >;
  
  // save memory, only preserve vertices_ and neighbours_
  struct facet_reduced
  {
    // each neighbouring facet lies against corresponding vertex and vice versa
    point_array vertices_; // D points (oriented)
    facet_array neighbours_; // D neighbouring facets
  };
  
  struct facet: facet_reduced // (d - 1)-dimensional face
  {
    point_list outside_; // if empty, then is convex hull's facet, else the first point (i.e. outside_.front()) is the furthest point from this facet
    point_forward_list coplanar_; // containing coplanar points and vertices of coplanar facets as well
    
    // equation of supporting hyperplane
    vector normal_; // components of normalized normal vector
    value_type D_; // distance from the origin to the hyperplane
    
    template< typename iterator >
    value_type
    distance(iterator const _point) const
    {
      using iterator_traits = std::iterator_traits< iterator >;
      static_assert(std::is_base_of< std::input_iterator_tag, typename iterator_traits::iterator_category >::value, "");
      return std::inner_product(normal_.cbegin(), normal_.cend(), _point, D_);
    }
    
  };
  
  // move information to facet base class, and drop facets that the origin is not at most second least idx
  void facet_make_reduce(const point_iterator &pcenter, unsigned long* nneighbor=nullptr)
  {
    std::unordered_set<long> nb_p_set; // count number of neighbors
    for(long source=facets_.size()-1; source>=0; --source)
    {
      // check if pcenter is smaller than the second largest vertex in facet
      int nsmaller = 0;
      facet & fc = facets_[source];
      for(auto &nb: fc.vertices_)
      {
        long rel_idx = nb-pcenter;
        nb_p_set.insert(rel_idx);
        if(rel_idx<0) ++nsmaller;
      }
      
      long lastfc = facets_.size()-1;
      if(nsmaller>1)
      {
        if(source<lastfc)
        {
          // remove this facet since pcenter is smaller than the second smallest vertex
          fc = std::move(facets_.back());
          for (small_uint_type const n : fc.neighbours_)
          {
            if(n<lastfc)
              replace_neighbour(n, (small_uint_type)lastfc, (small_uint_type)source);
          }
        }
        facets_.pop_back();
      }
      
    }
    if(nneighbor) *nneighbor=nb_p_set.size();
    while (!facets_.empty())
    {
      facets_reduced_.push_back({facets_.front().vertices_, facets_.front().neighbours_});
      facets_.pop_front();
    }
  }
  
  using facets = std::deque< facet >;
  using facets_reduced = std::deque< facet_reduced >;
  facets facets_;
  facets_reduced facets_reduced_;
  value_type
  cos_of_dihedral_angle(facet const & _first, facet const & _second) const
  {
    return std::inner_product(_first.normal_.cbegin(), _first.normal_.cend(), _second.normal_.cbegin(), zero);
  }
  
  private :
  
  void
  make_facet(facet & _facet,
             point_array const & _vertices,
             small_uint_type const _against,
             point_iterator const _apex,
             small_uint_type const _neighbour)
  {
    assert(_vertices.size() == D);
    _facet.vertices_ = _vertices;
    _facet.vertices_[_against] = _apex;
    //_facet.neighbours_.resize(D);
    _facet.neighbours_[_against] = _neighbour;
    //_facet.normal_.resize(D);
  }
  
  template< typename iterator >
  void
  make_facet(facet & _facet,
             iterator sbeg, // simplex
             small_uint_type const _vertex,
             bool const _swap)
  {
    using iterator_traits = std::iterator_traits< iterator >;
    static_assert(std::is_base_of< std::input_iterator_tag, typename iterator_traits::iterator_category >::value, "");
    static_assert(std::is_constructible< point_iterator, typename iterator_traits::value_type >::value, "");
    //_facet.vertices_.reserve(D);
    //_facet.neighbours_.reserve(D);
    for (small_uint_type v = 0, vb = 0; v <= D; ++v) {
      if (v != _vertex) {
        _facet.vertices_[vb] = *sbeg;
        _facet.neighbours_[vb] = v;
        ++vb;
      }
      ++sbeg;
    }
    if (_swap == (((D - _vertex) % 2) == 0)) {
      using std::swap;
      swap(_facet.vertices_.front(), _facet.vertices_.back());
      swap(_facet.neighbours_.front(), _facet.neighbours_.back());
    }
    //_facet.normal_.resize(D);
  }
  
  void
  reuse_facet(facet & _facet,
              point_array const & _vertices,
              small_uint_type const _against,
              point_iterator const _apex,
              small_uint_type const _neighbour)
  {
    //assert(_vertices.size() == D);
    _facet.vertices_ = _vertices;
    _facet.vertices_[_against] = _apex;
    //assert(_facet.neighbours_.size() == D);
    _facet.neighbours_[_against] = _neighbour;
    //assert(_facet.normal_.size() == D);
  }
  
  void
  copy_point(point_iterator const _from, vrow _to) const
  {
    std::copy_n((*_from).cbegin(), D, _to);
  }
  
  void
  subtract(vrow _minuend, crow _subtrahend) const
  {
    for (small_uint_type i = 0; i < D; ++i) {
      *_minuend++ -= *_subtrahend++;
    }
  }
  
  void
  gshift(vrow _augend, value_type const & _addend) const // shift Gaussian row
  {
    for (small_uint_type i = 0; i < D; ++i) {
      *_augend++ += _addend;
    }
    *_augend += _addend;
  }
  
  void
  divide(vrow _dividend, value_type const & _divisor) const
  {
    for (small_uint_type i = 0; i < D; ++i) {
      *_dividend++ /= _divisor;
    }
  }
  
  void
  subtract_and_assign(vrow _assignee, vrow _minuend, crow _subtrahend) const
  {
    for (small_uint_type i = 0; i < D; ++i) {
      *_assignee++ = (*_minuend++ -= *_subtrahend++);
    }
  }
  
  void
  multiply_and_add(vrow _assignee, crow _multiplicand, value_type const & _factor) const
  {
    for (small_uint_type i = 0; i < D; ++i) {
      *_assignee++ += (*_multiplicand++ * _factor);
    }
  }
  
  void
  scale_and_shift(vrow _multiplicand, crow _direction, value_type const & _factor) const
  {
    for (small_uint_type i = 0; i < D; ++i) {
      (*_multiplicand++ *= _factor) += *_direction++;
    }
  }
  
  void
  matrix_restore(small_uint_type const _identity)
  {
    for (small_uint_type c = 0; c < _identity; ++c) {
      std::copy_n(shadow_matrix_[c], D, matrix_[c]);
    }
    std::fill_n(matrix_[_identity], D, one);
    for (small_uint_type c = _identity + 1; c < D; ++c) {
      std::copy_n(shadow_matrix_[c], D, matrix_[c]);
    }
  }
  
  void
  matrix_restore()
  {
    for (small_uint_type c = 0; c < D; ++c) {
      std::copy_n(shadow_matrix_[c], D, matrix_[c]);
    }
  }
  
  void
  matrix_sqr(small_uint_type const _size)
  { // shadow_matrix_ = matrix_ * transposed matrix_
    assert(_size < D);
    for (small_uint_type r = 0; r < _size; ++r) {
      vrow const lhs_ = shadow_matrix_[r];
      crow const row_ = matrix_[r];
      for (small_uint_type c = 0; c < _size; ++c) {
        lhs_[c] = std::inner_product(row_, row_ + _size, matrix_[c], zero);
      }
    }
  }
  
  // based on LUP decomposition (complexity is (d^3 / 3 + d^2 / 2 - 5 * d / 6) vs (2 * d^3 / 3 + d^2 + d / 3 - 2) for QR decomposition via Householder reflections)
  value_type
  det(matrix const & _matrix, small_uint_type const _dimension)
  { // det_matrix_ contains lower unit triangular matrix and upper triangular at return
    assert(0 < _dimension);
    value_type det_ = one;
    std::copy_n(_matrix.cbegin(), _dimension, std::begin(det_matrix_));
    for (small_uint_type i = 0; i < _dimension; ++i) {
      vrow & mi_ = det_matrix_[i];
      small_uint_type pivot = i;
      {
        using std::abs;
        value_type max_ = abs(mi_[i]);
        small_uint_type j = i;
        while (++j < _dimension) {
          value_type y_ = abs(det_matrix_[j][i]);
          if (max_ < y_) {
            max_ = std::move(y_);
            pivot = j;
          }
        }
        if (!(eps < max_)) { // regular?
          det_ = zero; // singular
          break;
        }
      }
      if (pivot != i) {
        det_ = -det_; // each permutation flips sign of det
        std::swap(mi_, det_matrix_[pivot]);
      }
      value_type const & dia_ = mi_[i];
      det_ *= dia_; // det is multiple of diagonal elements
      small_uint_type j = i;
      while (++j < _dimension) {
        vrow const mj_ = det_matrix_[j];
        value_type & mji_ = mj_[i];
        mji_ /= dia_;
        small_uint_type k = i;
        while (++k < _dimension) {
          mj_[k] -= mji_ * mi_[k];
        }
      }
    }
    return det_;
  }
  
  // copy-paste from det(matrix const & _matrix, size_type const _dimension)
  // return det(det_matrix_, D);
  value_type gauss_eliminate(vector &x, vector &b)
  {
    value_type tmp, det_ = one;;
    small_int_type i, j;
    for (i = 0; i < D; ++i) {
      vrow & mi_ = det_matrix_[i];
      small_int_type pivot = i;
      {
        using std::abs;
        value_type max_ = abs(mi_[i]);
        j = i;
        while (++j < D) {
          value_type y_ = abs(det_matrix_[j][i]);
          if (max_ < y_) {
            max_ = std::move(y_);
            pivot = j;
          }
        }
      }
      if (pivot != i) {
        vrow & mj_ = det_matrix_[pivot];
        for(j=i; j<D; ++j)
        {
          tmp = mi_[j]; mi_[j] = mj_[j]; mj_[j] = tmp;
        }
        tmp = b[pivot]; b[pivot] = b[i]; b[i] = tmp;
        det_ = -det_;
      }
      value_type const & dia_ = mi_[i];
      det_ *= dia_;
      j = i;
      while (++j < D) {
        vrow const mj_ = det_matrix_[j];
        value_type & mji_ = mj_[i];
        mji_ /= dia_;
        small_int_type k = i;
        while (++k < D) {
          mj_[k] -= mji_ * mi_[k];
        }
        b[j] -= mji_ * b[i];
      }
    }
    while (--i >= 0) {
      vrow const mi_ = det_matrix_[i];
      value_type &bi = b[i];
      j = D;
      while (--j > i)
        bi -= x[j] * mi_[j];
      x[i] = bi / mi_[i];
    }
    return det_;
  }
  
  void
  set_hyperplane_equation(facet & _facet)
  {
    vector x, b;
    for(small_uint_type i=0; i<D; ++i)
    {
      std::copy_n(_facet.vertices_[i]->cbegin(), D, det_matrix_[i]);
      b[i] = one;
    }
    value_type det_ = gauss_eliminate(x, b);
    value_type N = zero;
    for (small_uint_type i = 0; i < D; ++i) {
      N += x[i] * x[i];
    }
    N = 1.0/std::sqrt(std::move(N));
    if(0. < det_)
    {
      for (small_uint_type i = 0; i < D; ++i) {
       _facet.normal_[i] = x[i]*N;
      }
      _facet.D_ = -N;
    }
    else
    {
      for (small_uint_type i = 0; i < D; ++i) {
       _facet.normal_[i] = -x[i]*N;
      }
      _facet.D_ = N;
    }
    assert(_facet.distance(inner_point_) < zero);
  }
  
  bool
  orthonormalize(point_list const & _affine_space,
                 small_uint_type const _rank,
                 crow const _origin)
  {
    assert(!(D < _rank));
    assert(!(_affine_space.size() < _rank));
    auto vertex = std::begin(_affine_space);
    for (small_uint_type r = 0; r < _rank; ++r) { // affine space -> vector space
      vrow const row_ = shadow_matrix_[r];
      copy_point(*vertex, row_);
      subtract(row_, _origin);
      ++vertex;
    }
    for (small_uint_type i = 0; i < _rank; ++i) { // Householder transformation
      value_type sum_ = zero;
      vrow const qri_ = shadow_matrix_[i];
      for (small_uint_type k = i; k < D; ++k) {
        value_type const & qrik_ = qri_[k];
        sum_ += qrik_ * qrik_;
      }
      using std::sqrt;
      value_type norm_ = sqrt(sum_);
      if (!(eps < norm_)) {
        return false;
      }
      value_type & qrii_ = qri_[i];
      if (qrii_ < zero) {
        norm_ = -norm_;
      }
      value_type factor_ = sqrt(std::move(sum_) + qrii_ * norm_);
      if (!(eps < factor_)) {
        return false;
      }
      qrii_ += std::move(norm_);
      for (small_uint_type k = i; k < D; ++k) {
        qri_[k] /= factor_;
      }
      small_uint_type j = i;
      while (++j < _rank) {
        vrow const qrj_ = shadow_matrix_[j];
        value_type s_ = zero;
        for (small_uint_type k = i; k < D; ++k) {
          s_ += qri_[k] * qrj_[k];
        }
        for (small_uint_type k = i; k < D; ++k) {
          qrj_[k] -= qri_[k] * s_;
        }
      }
    } // shadow_matrix_ is packed QR
    return true;
  }
  
  void
  forward_transformation(small_uint_type const _rank) // calculation of Q
  {
    assert(!(D < _rank));
    for (small_uint_type i = 0; i < _rank; ++i) {
      vrow const qi_ = matrix_[i];
      std::fill_n(qi_, D, zero);
      qi_[i] = one;
      small_uint_type j = _rank;
      while (0 < j) {
        --j;
        vrow const qrj_ = shadow_matrix_[j]; // containing packed QR
        value_type s_ = zero;
        for (small_uint_type k = j; k < D; ++k) {
          s_ += qrj_[k] * qi_[k];
        }
        for (small_uint_type k = j; k < D; ++k) {
          qi_[k] -= qrj_[k] * s_;
        }
      }
    } // matrix_ is Q
  }
  
  bool
  steal_best(point_list & _basis)
  { // set moves a point which is furthest from affine subspace formed by points of "_basis" set from "outside_" set to "_basis"
    assert(!_basis.empty());
    small_uint_type const rank_ = (small_uint_type)_basis.size() - 1;
    assert(rank_ < D);
    vrow const origin_ = matrix_[rank_];
    copy_point(_basis.back(), origin_);
    if (!orthonormalize(_basis, rank_, origin_)) {
      return false;
    }
    forward_transformation(rank_);
    vrow const projection_ = shadow_matrix_.back();
    vrow const apex_ = shadow_matrix_.front();
    value_type distance_ = zero; // square of distance to the subspace
    auto oend = outside_.end();
    auto furthest = oend;
    for (auto it = outside_.begin(); it != oend; ++it) {
      copy_point(*it, apex_);
      subtract_and_assign(projection_, apex_, origin_); // turn translated space into vector space then project onto orthogonal subspace
      for (small_uint_type i = 0; i < rank_; ++i) {
        crow const qi_ = matrix_[i];
        multiply_and_add(projection_, qi_, -std::inner_product(qi_, qi_ + D, apex_, zero));
      }
      value_type d_ = std::inner_product(projection_, projection_ + D, projection_, zero);
      if (distance_ < d_) {
        distance_ = std::move(d_);
        furthest = it;
      }
    }
    if (furthest == oend) {
      return false;
    }
    _basis.splice(_basis.end(), outside_, furthest);
    return true;
  }
  
  facet_container removed_facets_;
  
  std::pair< facet &, small_uint_type const >
  add_facet(point_array const & _vertices,
            small_uint_type const _against,
            point_iterator const _apex,
            small_uint_type const _neighbour)
  {
    if (removed_facets_.empty()) {
      small_uint_type const f = (small_uint_type)facets_.size();
      facets_.emplace_back();
      facet & facet_ = facets_.back();
      make_facet(facet_, _vertices, _against, _apex, _neighbour);
      return {facet_, f};
    } else {
      small_uint_type const f = removed_facets_.back();
      removed_facets_.pop_back();
      facet & facet_ = facets_[f];
      reuse_facet(facet_, _vertices, _against, _apex, _neighbour);
      return {facet_, f};
    }
  }
  
  using ranking = std::multimap< value_type, small_uint_type >;
  using ranking_meta = std::unordered_map< small_uint_type, typename ranking::iterator >;
  
  ranking ranking_;
  ranking_meta ranking_meta_;
  
  void
  rank(value_type && _orientation,
       small_uint_type const f)
  {
    if (eps < _orientation) {
      ranking_meta_.emplace(f, ranking_.emplace(std::move(_orientation), f));
    }
  }
  
  void
  unrank(small_uint_type const f)
  {
    auto const r = ranking_meta_.find(f);
    if (r != std::end(ranking_meta_)) {
      ranking_.erase(r->second);
      ranking_meta_.erase(r);
    }
    removed_facets_.push_back(f);
  }
  
  point_list outside_;
  
  value_type
  partition(facet & _facet)
  {
    value_type distance_ = zero;
    auto it = outside_.begin();
    auto oend = outside_.end();
    while (it != oend) {
      auto const next = std::next(it);
      value_type d_ = _facet.distance((**it).cbegin());
      if (eps < d_) {
        if (distance_ < d_) {
          distance_ = std::move(d_);
          _facet.outside_.splice(_facet.outside_.begin(), outside_, it);
        } else {
          _facet.outside_.splice(_facet.outside_.end(), outside_, it);
        }
      } else if (!(d_ < -eps)) {
        _facet.coplanar_.push_front(*it);
      }
      it = next;
    }
    return distance_;
  }
  
  small_uint_type
  get_best_facet() const
  {
    assert(ranking_meta_.size() == ranking_.size());
    return std::prev(ranking_.cend())->second;
  }
  
  void
  replace_neighbour(small_uint_type const f,
                    small_uint_type const _from,
                    small_uint_type const _to)
  {
    if (_from != _to) {
      for (small_uint_type & n : facets_[f].neighbours_) {
        if (n == _from) {
          n = _to;
          return;
        }
      }
    }
  }

  struct ridge
  {
    
    facet & facet_;
    small_uint_type const f;
    small_uint_type const v;
    size_type const hash_;
    
    bool
    operator == (ridge const & _rhs) const noexcept
    {
      point_iterator const lskip = facet_.vertices_[v];
      point_iterator const rskip = _rhs.facet_.vertices_[_rhs.v];
      for (point_iterator const & l : facet_.vertices_) {
        if (l != lskip) {
          bool found_ = false;
          for (point_iterator const & r : _rhs.facet_.vertices_) {
            if (r != rskip) {
              if (l == r) {
                found_ = true; // O(D^2) expensive
                break;
              }
            }
          }
          if (!found_) {
            return false;
          }
        }
      }
      return true;
    }
    
  };
  
  struct ridge_hash
  {
    
    size_type
    operator () (ridge const & _ridge) const noexcept
    {
      return _ridge.hash_;
    }
    
  };
  
  std::unordered_set< ridge, ridge_hash > unique_ridges_;
  std::hash< typename std::iterator_traits< point_iterator >::value_type const * > point_hash_;
  std::array< size_type, D > vertices_hashes_;
  
  void
  find_adjacent_facets(facet & _facet,
                       small_uint_type const f,
                       small_uint_type const _skip)
  {
    size_type ridge_hash_ = 0;
    for (small_uint_type v = 0; v < D; ++v) {
      if (v != _skip) {
        ridge_hash_ ^= (vertices_hashes_[v] = point_hash_(std::addressof(*_facet.vertices_[v])));
      }
    }
    for (small_uint_type v = 0; v < D; ++v) {
      if (v != _skip) { // neighbouring facet against apex (_skip-indexed) is known atm
        auto const position = unique_ridges_.insert({_facet, f, v, (ridge_hash_ ^ vertices_hashes_[v])});
        if (!position.second) {
          ridge const & ridge_ = *position.first;
          ridge_.facet_.neighbours_[ridge_.v] = f;
          _facet.neighbours_[v] = ridge_.f;
          unique_ridges_.erase(position.first);
        }
      }
    }
  }
  
  using facet_unordered_set = std::unordered_set< small_uint_type >;
  
  facet_unordered_set visited_;
  facet_unordered_set visible_;
  
  // function process_visibles is now recursion-free
  struct visible_status
  {
    small_uint_type f_;
    small_uint_type v_;
    visible_status(small_uint_type f): f_(f), v_(0) {}
  };
  
  // check if the distance of point _apex to facet facets_[f0] is positive
  inline bool is_point_in_facet(small_uint_type const f,
                                point_iterator const _apex,
                                bool &last_flag)
  {
    assert(!(visited_.size() < visible_.size()));
    if (!visited_.insert(f).second) {
      last_flag = (visible_.count(f) != 0);
      return false;
    }
    facet & facet_ = facets_[f];
    if (!(zero < facet_.distance((*_apex).cbegin()))) {
      last_flag = false;
      return false;
    }
    visible_.insert(f);
    outside_.splice(outside_.end(), std::move(facet_.outside_));
    facet_.coplanar_.clear();
    return true;
  }
  
  bool
  process_visibles(facet_container & _newfacets,
                   small_uint_type const f0,
                   point_iterator const _apex) // traverse the graph of visible facets
  {
    std::stack<visible_status> fs_;
    bool last_flag = true;
    if(is_point_in_facet(f0, _apex, last_flag)) fs_.push(visible_status(f0));
    while(!fs_.empty())
    {
      auto &pr = fs_.top();
      small_uint_type f = pr.f_;
      facet & facet_ = facets_[f];
      small_uint_type v = pr.v_; // v = [0, D]
      if(v < D)
      {
        small_uint_type n = facet_.neighbours_[v];
        if(is_point_in_facet(n, _apex, last_flag))
        {
          fs_.push(visible_status(n));
        }
        else
        {
          if(!last_flag)
          {
            auto const newfacet = add_facet(facet_.vertices_, v, _apex, n);
            set_hyperplane_equation(newfacet.first);
            _newfacets.push_back(newfacet.second);
            replace_neighbour(n, f, newfacet.second);
            find_adjacent_facets(newfacet.first, newfacet.second, v);
          }
          ++pr.v_;
        }
      }
      else // == D
      {
        unrank(f);
        last_flag = true;
        fs_.pop();
      }
    }
    return last_flag;
  }
  
  
  void
  compactify()
  {
    small_uint_type source = (small_uint_type)facets_.size();
    assert(removed_facets_.size() < source);
    assert(D < source - removed_facets_.size());
    assert(ranking_.size() == ranking_meta_.size());
    assert(!(source < ranking_.size()));
    auto const rend = std::end(ranking_meta_);
    std::sort(removed_facets_.rbegin(), removed_facets_.rend());
    for (small_uint_type const destination : removed_facets_) {
      assert(!(source < destination));
      if (destination != --source) {
        facet & facet_ = facets_[destination];
        facet_ = std::move(facets_.back());
        for (small_uint_type const n : facet_.neighbours_) {
          replace_neighbour(n, source, destination);
        }
        auto const r = ranking_meta_.find(source);
        if (r != rend) {
          r->second->second = destination;
          ranking_meta_.emplace(destination, std::move(r->second));
          ranking_meta_.erase(r);
        }
      }
      facets_.pop_back();
    }
    removed_facets_.clear();
  }
  
  bool
  check_local_convexity(facet const & facet_,
                        small_uint_type const f) const
  {
    assert(&facets_[f] == &facet_);
    for (small_uint_type const n : facet_.neighbours_) {
      facet const & neighbour_ = facets_[n];
      if (cos_of_dihedral_angle(facet_, neighbour_) < one) { // avoid roundoff error
        for (small_uint_type v = 0; v < D; ++v) {
          if (neighbour_.neighbours_[v] == f) { // vertex v of neigbour_ facet is opposite to facet_
            value_type const distance_ = facet_.distance((*neighbour_.vertices_[v]).cbegin());
            if (eps < distance_) {
              return false; // facet is not locally convex at ridge, common for facet_ and neighbour_ facets
            } else {
              break;
            }
          }
        }
      }
    }
    return true;
  }
  
  public :
  
  template< typename iterator >
  value_type
  hypervolume(iterator first,
              iterator const last) // hypervolume of parallelotope spanned on vectors from last vertex (vlast) to all the vertices lies in [vfirst, vlast)
  {
    using iterator_traits = std::iterator_traits< iterator >;
    static_assert(std::is_base_of< std::input_iterator_tag, typename iterator_traits::iterator_category >::value, "");
    static_assert(std::is_constructible< point_iterator, typename iterator_traits::value_type >::value, "");
    if (first == last) {
      return zero;
    }
    vrow const origin_ = shadow_matrix_.back();
    copy_point(*last, origin_);
    small_uint_type rank_ = 0;
    while (first != last) { // affine space -> vector space
      assert(rank_ < D);
      vrow const row_ = matrix_[rank_];
      copy_point(*first, row_);
      subtract(row_, origin_);
      ++rank_;
      ++first;
    }
    if (rank_ == D) {
      return det(matrix_, D); // oriented hypervolume
    } else {
      matrix_sqr(rank_);
      using std::sqrt;
      return sqrt(det(shadow_matrix_, rank_)); // non-oriented rank_-dimensional measure
    }
  }
  
  void
  add_points(point_iterator beg,
             point_iterator const end) // [beg; end)
  {
    while (beg != end) {
      outside_.push_back(beg);
      ++beg;
    }
  }
  
  template< typename iterator >
  void
  add_points(iterator const beg,
             iterator const end) // [beg; end)
  {
    using iterator_traits = std::iterator_traits< iterator >;
    static_assert(std::is_base_of< std::input_iterator_tag, typename iterator_traits::iterator_category >::value, "");
    static_assert(std::is_constructible< point_iterator, typename iterator_traits::value_type >::value, "");
    std::copy(beg, end, std::back_inserter(outside_));
  }
  
  point_list
  get_affine_basis()
  {
    assert(facets_.empty());
    point_list basis_;
    basis_.splice(basis_.end(), outside_, std::begin(outside_));
    if (!steal_best(basis_)) {
      return basis_; // can't find affinely independent second point
    }
    outside_.splice(outside_.begin(), basis_, std::begin(basis_)); // reject first point to rejudge it
    for (small_uint_type i = 0; i < D; ++i) {
      if (!steal_best(basis_)) {
        return basis_; // can't find (i + 2) affinely independent points
      }
    }
    return basis_;
  }
  
  template< typename iterator >
  value_type
  create_initial_simplex(iterator const first,
                         iterator const last) // [bfirst; blast]
  {
    using iterator_traits = std::iterator_traits< iterator >;
    static_assert(std::is_base_of< std::forward_iterator_tag, typename iterator_traits::iterator_category >::value, "");
    static_assert(std::is_constructible< point_iterator, typename iterator_traits::value_type >::value, "");
    assert(static_cast< small_uint_type >(std::distance(first, last)) == D);
    assert(facets_.empty());
    {
      copy_point(*last, inner_point_);
      auto it = first;
      while (it != last) {
        auto x = (**it).cbegin();
        for (small_uint_type i = 0; i < D; ++i) {
          inner_point_[i] += *x;
          ++x;
        }
        ++it;
      }
      divide(inner_point_, value_type(D + 1));
    }
    value_type const volume_ = hypervolume(first, last);
    bool const swap_ = (volume_ < zero);
    for (small_uint_type f = 0; f <= D; ++f) {
      facets_.emplace_back();
      facet & facet_ = facets_.back();
      make_facet(facet_, first, f, swap_);
      set_hyperplane_equation(facet_);
      rank(partition(facet_), f);
    }
    outside_.clear();
    assert(check());
    return volume_;
  }
  
  // Barber, C. B., D.P. Dobkin, and H.T. Huhdanpaa, 1995. "The Quickhull Algorithm for Convex Hulls", ACM Transactions on Mathematical Software.
  void
  create_convex_hull()
  {
    assert(facets_.size() == D + 1);
    assert(removed_facets_.empty());
    facet_container newfacets_;
    while (!ranking_.empty()) {
      small_uint_type const f = get_best_facet();
      point_list & o_ = facets_[f].outside_;
      assert(!o_.empty());
      point_iterator const apex = std::move(o_.front());
      o_.pop_front();
      if (!process_visibles(newfacets_, f, apex)) {
        assert(false);
      }
      visited_.clear();
      visible_.clear();
      assert(unique_ridges_.empty());
      for (small_uint_type const n : newfacets_) {
        facet & facet_ = facets_[n];
        // Time consumming check. Invoke when necessary
        // assert(check_local_convexity(facet_, n));
        rank(partition(facet_), n);
      }
      newfacets_.clear();
      outside_.clear();
      //assert((compactify(), check()));
    }
    assert(ranking_meta_.empty());
    compactify();
  }
  
  // Kurt Mehlhorn, Stefan Näher, Thomas Schilz, Stefan Schirra, Michael Seel, Raimund Seidel, and Christian Uhrig.
  // Checking geometric programs or verification of geometric structures. In Proc. 12th Annu. ACM Sympos. Comput. Geom., pages 159–165, 1996.
  bool
  check() const
  {
    assert(D < facets_.size());
    small_uint_type facets_count_ = 0;
    for (facet const & facet_ : facets_) { // check local convexity of all the facets
      if (!check_local_convexity(facet_, facets_count_)) {
        return false;
      }
      ++facets_count_;
    }
    facet const & first_ = facets_.front();
    {
      value_type const distance_ = first_.distance(inner_point_);
      if (!(distance_ < zero)) {
        return false; // inner point is not on negative side of the first facet, therefore structure is not convex
      }
    }
    std::vector<value_type> memory_(D * (4 + D), zero);
    vrow centroid_ = memory_.data();
    vrow const ray_ = centroid_;
    centroid_ += D;
    for (point_iterator const & v : first_.vertices_) {
      auto x = (*v).cbegin();
      for (small_uint_type i = 0; i < D; ++i) {
        ray_[i] += *x;
        ++x;
      }
    }
    divide(ray_, value_type(D));
    subtract(ray_, inner_point_);
    {
      value_type const dot_product_ = std::inner_product(ray_, ray_ + D, first_.normal_.data(), zero);
      if (!(zero < dot_product_)) { // ray is parallel to the plane or directed away from the plane
        return false;
      }
    }
    matrix g_; // storage (d * (d + 1)) for Gaussian elimination with partial pivoting
    for (vrow & row_ : g_) {
      row_ = centroid_;
      centroid_ += (D + 1);
    }
    vrow const intersection_point_ = centroid_;
    centroid_ += D;
    assert(centroid_ + D == &memory_.back() + 1);
    for (small_uint_type f = 1; f < facets_count_; ++f) {
      using std::abs;
      facet const & facet_ = facets_[f];
      value_type const numerator_ = facet_.distance(inner_point_);
      if (!(numerator_ < zero)) {
        return false; // inner point is not on negative side of all the facets, i.e. structure is not convex
      }
      value_type const denominator_ = std::inner_product(ray_, ray_ + D, facet_.normal_.data(), zero);
      if (!(zero < denominator_)) { // ray is parallel to the plane or directed away from the plane
        continue;
      }
      std::copy_n(ray_, D, intersection_point_);
      scale_and_shift(intersection_point_, inner_point_, -(numerator_ / denominator_));
      for (small_uint_type v = 0; v < D; ++v) {
        auto beg = (*facet_.vertices_[v]).cbegin();
        for (small_uint_type r = 0; r < D; ++r) {
          g_[r][v] = *beg;
          ++beg;
        }
      }
      for (small_uint_type r = 0; r < D; ++r) {
        vrow const gr_ = g_[r];
        centroid_[r] = -std::accumulate(gr_, gr_ + D, zero) / value_type(D);
        gr_[D] = intersection_point_[r];
      }
      for (small_uint_type r = 0; r < D; ++r) {
        vrow const gr_ = g_[r];
        value_type & x_ = centroid_[r];
        gshift(gr_, x_);
        //assert(!(eps * value_type(D) < std::accumulate(gr_, gr_ + D, zero))); // now center of the facet coincides with the origin, but no one vertex does
        auto const bounding_box = std::minmax_element(gr_, gr_ + D);
        x_ = *bounding_box.second - *bounding_box.first;
        if (!(eps * value_type(D) < x_)) {
          x_ = one;
        }
      }
      for (small_uint_type r = 0; r < D; ++r) {
        vrow const gr_ = g_[r];
        gshift(gr_, centroid_[r]);
      }
      for (small_uint_type i = 0; i < D; ++i) { // Gaussian elimination
        vrow & gi_ = g_[i];
        value_type max_ = abs(gi_[i]);
        small_uint_type pivot = i;
        {
          small_uint_type p = i;
          while (++p < D) {
            value_type y_ = abs(g_[p][i]);
            if (max_ < y_) {
              max_ = std::move(y_);
              pivot = p;
            }
          }
        }
        assert(eps < max_); // vertex must not match the origin after above transformations
        if (pivot != i) {
          std::swap(gi_, g_[pivot]);
        }
        value_type & gii_ = gi_[i];
        for (small_uint_type j = i + 1; j < D; ++j) {
          vrow const gj_ = g_[j];
          value_type & gji_ = gj_[i];
          gji_ /= gii_;
          for (small_uint_type k = i + 1; k <= D; ++k) {
            gj_[k] -= gji_ * gi_[k];
          }
          gji_ = zero;
        }
      } // g_ is upper triangular now
      bool in_range_ = true;
      {
        small_uint_type i = D;
        while (0 < i) {
          --i;
          vrow const gi_ = g_[i];
          value_type & xi_ = gi_[D];
          for (small_uint_type j = i + 1; j < D; ++j) {
            xi_ -= gi_[j] * g_[j][D];
          }
          value_type const & gii_ = gi_[i];
          assert(eps < abs(gii_)); // vertex must not match the origin
          xi_ /= gii_;
          if ((xi_ < zero) || (one < xi_)) {
            in_range_ = false; // barycentric coordinate does not lie in [0;1] interval => miss
            break;
          }
        }
      }
      if (in_range_) {
        return false; // hit
      }
    }
    return true;
  }
// debug functions
public:
  void checksizes()
  {
    std::cout << "quick_hull_size: " << sizeof(*this) << std::endl;
    std::cout << "facet_reduced: " << sizeof(facet_reduced) << std::endl;
    std::cout << "facet: " << sizeof(facet) << std::endl;
    std::cout << "  vertices: " << sizeof(facet::vertices_) << std::endl;
    std::cout << "  neighbors: " << sizeof(facet::neighbours_) << std::endl;
    std::cout << "  outside: " << sizeof(facet::outside_) << std::endl;
    std::cout << "  coplanar: " << sizeof(facet::coplanar_) << std::endl;
    std::cout << "  normal: " << sizeof(facet::normal_) << std::endl;
    std::cout << "  D_: " << sizeof(facet::D_) << std::endl;
    std::cout << "ridge: " << sizeof(ridge) << std::endl;
    std::cout << "facet_container: " << sizeof(removed_facets_) << std::endl;
    std::cout << "visible_status: " << sizeof(visible_status) << std::endl;
    std::cout << "point_iterator: " << sizeof(point_iterator) << std::endl;
  }
};
