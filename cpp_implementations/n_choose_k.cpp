#include <stdio.h>
#include <stddef.h>
#include <stdint.h>
#include <iosfwd>
#include <cstddef>
#include <algorithm>


/* Calculates number of combinations of n things taken k at a time */
/* The code has been retrieved from https://stackoverflow.com/a/4776276 */
unsigned long long n_choose_k(const unsigned long long& n,
    const unsigned long long& k)
{
    if (n < k) return 0;
    if (0 == n) return 0;
    if (0 == k) return 1;
    if (n == k) return 1;
    if (1 == k) return n;
    typedef unsigned long long value_type;
    value_type* table = new value_type[static_cast<std::size_t>(n * n)];
    std::fill_n(table, n * n, 0);
    class n_choose_k_impl
    {
    public:

        n_choose_k_impl(value_type* table, const value_type& dimension)
            : table_(table),
            dimension_(dimension)
        {}

        inline value_type& lookup(const value_type& n, const value_type& k)
        {
            return table_[dimension_ * n + k];
        }

        inline value_type compute(const value_type& n, const value_type& k)
        {
            if ((0 == k) || (k == n))
                return 1;
            value_type v1 = lookup(n - 1, k - 1);
            if (0 == v1)
                v1 = lookup(n - 1, k - 1) = compute(n - 1, k - 1);
            value_type v2 = lookup(n - 1, k);
            if (0 == v2)
                v2 = lookup(n - 1, k) = compute(n - 1, k);
            return v1 + v2;
        }

        value_type* table_;
        value_type dimension_;
    };
    value_type result = n_choose_k_impl(table, n).compute(n, k);
    delete[] table;
    return result;
}