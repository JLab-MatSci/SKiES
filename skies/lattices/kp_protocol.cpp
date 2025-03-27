/*-----------------------------------------------------------------------
    * SKiES - Solver of Kinetic Equation for Solids
    * 
    * (C) 2025 Galtsov Ilya, Fokin Vladimir, Minakov Dmitry, Levashov Pavel (JIHT RAS)
    *
    * SKiES may only be utilized for non-profit research.
    * Citing appropriate sources is required when using SKiES.
    * 
    * Distribution of this file is permitted by the GNU General Public License.
    * Examine the `LICENSE' file located in the current distribution's root directory.
------------------------------------------------------------------------- */
#include <cmath>
#include <cassert>
#include <algorithm>
#include <float.h>
#include <iostream>
#include <numeric>

#include <skies/common/alg.h>
#include <skies/lattices/kp_protocol.h>

using namespace skies::arrays;

namespace skies {

bool operator== (const KPprotocol& kp1, const KPprotocol& kp2)
{
    return (kp1.n1_ == kp2.n1_) &&
           (kp1.n2_ == kp2.n2_) &&
           (kp1.n3_ == kp2.n3_) &&
           (kp1.nkpt_ == kp2.nkpt_) &&
           (kp1.grid_ == kp2.grid_) &&
           (kp1.igrid_ == kp2.igrid_) &&
           (kp1.kprange_ == kp2.kprange_) &&
           (kp1.vd_from_ind_ == kp2.vd_from_ind_) &&
           (kp1.vi_from_ind_ == kp2.vi_from_ind_);
}

KPprotocol::KPprotocol() {}

KPprotocol::KPprotocol(size_t n1, size_t n2, size_t n3)
    : n1_(n1), n2_(n2), n3_(n3)
    , nkpt_(n1 * n2 * n3)
    , grid_(nkpt_, array1D(3, 0.0))
    , igrid_(nkpt_, std::vector<int>(3, 0))
{
    std::vector<int> vi1(n1);
    std::vector<int> vi2(n2);
    std::vector<int> vi3(n3);
    std::iota(vi1.begin(), vi1.end(), 1);
    std::iota(vi2.begin(), vi2.end(), 1);
    std::iota(vi3.begin(), vi3.end(), 1);

    // translate to create a zero-symmetric interval, step = 1
    std::for_each(vi1.begin(), vi1.end(), [n1] (int& i) { i = 2 * i - n1 - 1; return i; });
    std::for_each(vi2.begin(), vi2.end(), [n2] (int& i) { i = 2 * i - n2 - 1; return i; });
    std::for_each(vi3.begin(), vi3.end(), [n3] (int& i) { i = 2 * i - n3 - 1; return i; });

    // MP-scheme
    int cnt{ 0 };
    for (auto i3 : vi3)
        for (auto i1 : vi1)
            for (auto i2 : vi2)
            {
                std::vector<int>& ki = igrid_[cnt];
                std::vector<double>& k = grid_[cnt];
                ki[0] = i1; k[0] = 0.5 * static_cast<double>(i1) / n1;
                ki[1] = i2; k[1] = 0.5 * static_cast<double>(i2) / n2;
                ki[2] = i3; k[2] = 0.5 * static_cast<double>(i3) / n3;
                // fill in map
                auto N = find_ind_from_vi(ki);
                kprange_.push_back(N);
                vi_from_ind_[N] = ki;
                vd_from_ind_[N] = k;
                cnt++;
            }
}

std::tuple<size_t, size_t, size_t> KPprotocol::mesh() const
{
    return std::make_tuple(n1_, n2_, n3_);
}

std::vector<size_t> KPprotocol::local_subcell(size_t ind) const
{
    auto ki = vi_from_ind_.at(ind);
    std::vector<size_t> subcell_i;
    std::vector<std::vector<int>> subcell_v(8);

    // indices correspond to ones preseneted on fig. 5 of Bloechl, Jepsen and Andersen prb 49.23 (1994): 16223.
    subcell_v[0] = ki;
    subcell_v[1] = ki + std::vector<int>{2, 0, 0};
    subcell_v[2] = ki + std::vector<int>{0, 2, 0};
    subcell_v[3] = ki + std::vector<int>{2, 2, 0};
    subcell_v[4] = ki + std::vector<int>{0, 0, 2};
    subcell_v[5] = ki + std::vector<int>{2, 0, 2};
    subcell_v[6] = ki + std::vector<int>{0, 2, 2};
    subcell_v[7] = ki + std::vector<int>{2, 2, 2};

    auto max = std::max_element(ki.begin(), ki.end());
    assert(max != ki.end());
    if (*max + 1 == static_cast<int>(n1_) || *max + 1 == static_cast<int>(n2_) || *max + 1 == static_cast<int>(n3_))
        this->apply_pbc(subcell_v);
    
    for (auto&& k : subcell_v)
        subcell_i.push_back(this->find_ind_from_vi(k));

    return subcell_i;
}

void KPprotocol::apply_pbc(std::vector<std::vector<int>>& subcell_v) const
{
    for (auto&& k : subcell_v)
    {
        if (k[0] > static_cast<int>(n1_ - 1))
            k[0] -= 2 * static_cast<int>(n1_);
        if (k[1] > static_cast<int>(n2_ - 1))
            k[1] -= 2 * static_cast<int>(n2_);
        if (k[2] > static_cast<int>(n3_ - 1))
            k[2] -= 2 * static_cast<int>(n3_);
    }
}

size_t KPprotocol::find_ind_from_vi(const std::vector<int>& vi) const
{
    auto i = (vi[0] + n1_ - 1) / 2;
    auto j = (vi[1] + n2_ - 1) / 2;
    auto k = (vi[2] + n3_ - 1) / 2;
    return n2_ * i + j + n1_ * n2_ * k;
}

std::vector<int> KPprotocol::find_vi_from_ind(size_t ind) const
{
    return vi_from_ind_.at(ind);
}

array1D KPprotocol::find_vd_from_ind(size_t ind) const
{
    return vd_from_ind_.at(ind);
}

const std::vector<size_t>& KPprotocol::range() const &
{
    return kprange_;
}

} // skies
