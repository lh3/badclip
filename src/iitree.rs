//! Implicit interval tree (a Rust port of `IITree.h`, from cgranges).
//!
//! Suppose there are `N = 2^(K+1) - 1` sorted numbers in an array `a[]`. They
//! implicitly form a complete binary tree of height `K+1`. We consider leaves
//! to be at level 0. The binary tree has the following properties:
//!
//! 1. The lowest `k-1` bits of nodes at level `k` are all 1. The `k`-th bit is
//!    0. The first node at level `k` is indexed by `2^k - 1`. The root of the
//!    tree is indexed by `2^K - 1`.
//!
//! 2. For a node `x` at level `k`, its left child is `x - 2^(k-1)` and the
//!    right child is `x + 2^(k-1)`.
//!
//! 3. For a node `x` at level `k`, it is a left child if its `(k+1)`-th bit is
//!    0. Its parent node is `x + 2^k`. Similarly, if the `(k+1)`-th bit is 1,
//!    `x` is a right child and its parent is `x - 2^k`.
//!
//! 4. For a node `x` at level `k`, there are `2^(k+1) - 1` nodes in the
//!    subtree descending from `x`, including `x`. The left-most leaf is
//!    `x & !(2^k - 1)` (masking the lowest `k` bits to 0).
//!
//! When numbers can't fill a complete binary tree, the parent of a node may not
//! be present in the array. The implementation here still mimics a complete
//! tree, though getting the special casing right is a little complex.
//!
//! As a sorted array can be considered as a binary search tree, we can
//! implement an interval tree on top of the idea. We only need to record, for
//! each node, the maximum end in the subtree descending from the node.
//!
//! Intervals are half-open `[st, en)`; a query `[st, en)` overlaps interval
//! `i` iff `a[i].st < en && st < a[i].en`. `S` is the coordinate scalar type
//! (any `Copy + PartialOrd`, e.g. `i64`/`u32`), `T` is the payload attached to
//! each interval.
//!
//! ```ignore
//! let mut t = IITree::new();
//! t.add(10, 20, "a");
//! t.add(15, 30, "b");
//! t.index();
//! for i in t.overlap(18, 19) {
//!     println!("{}-{} {}", t.start(i), t.end(i), t.data(i));
//! }
//! ```

use std::cmp::Ordering;

struct Interval<S, T> {
    st: S,
    en: S,
    max: S,
    data: T,
}

/// An implicit interval tree over a sorted `Vec` of intervals.
pub struct IITree<S, T> {
    a: Vec<Interval<S, T>>,
    /// Level of the root, `None` for an empty tree.
    max_level: Option<usize>,
    /// Set by [`IITree::index`]; cleared by [`IITree::add`].
    indexed: bool,
}

impl<S: Copy + PartialOrd, T> Default for IITree<S, T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<S: Copy + PartialOrd, T> IITree<S, T> {
    /// Creates an empty tree. Call [`IITree::add`] for every interval, then
    /// [`IITree::index`] once before querying.
    pub fn new() -> Self {
        IITree {
            a: Vec::new(),
            max_level: None,
            indexed: false,
        }
    }

    /// Adds the interval `[st, en)` with payload `data`. Invalidates the
    /// index; [`IITree::index`] must be called again before querying.
    pub fn add(&mut self, st: S, en: S, data: T) {
        self.a.push(Interval {
            st,
            en,
            max: en,
            data,
        });
        self.indexed = false;
    }

    /// Sorts the intervals by start and computes the per-node subtree maxima.
    /// After indexing, interval indices `0..len()` are in start order.
    pub fn index(&mut self) {
        // Stable sort so that ties keep insertion order.
        self.a
            .sort_by(|x, y| x.st.partial_cmp(&y.st).unwrap_or(Ordering::Equal));
        self.max_level = Self::index_core(&mut self.a);
        self.indexed = true;
    }

    /// Returns the level of the root, or `None` if `a` is empty.
    fn index_core(a: &mut [Interval<S, T>]) -> Option<usize> {
        let n = a.len();
        if n == 0 {
            return None;
        }
        // last_i points to the rightmost node in the tree; last is its max.
        let mut last_i = 0usize;
        let mut last = a[0].en;
        // Leaves (i.e. at level 0).
        let mut i = 0;
        while i < n {
            a[i].max = a[i].en;
            last_i = i;
            last = a[i].en;
            i += 2;
        }
        // Process internal nodes in the bottom-up order.
        let mut k = 1usize;
        while (1usize << k) <= n {
            let x = 1usize << (k - 1);
            let i0 = (x << 1) - 1; // the first node at level k
            let step = x << 2;
            let mut i = i0;
            while i < n {
                // Traverse all nodes at level k.
                let el = a[i - x].max; // max value of the left child
                let er = if i + x < n { a[i + x].max } else { last }; // of the right child
                let mut e = a[i].en;
                if el > e {
                    e = el;
                }
                if er > e {
                    e = er;
                }
                a[i].max = e; // set the max value for node i
                i += step;
            }
            // last_i now points to the parent of the original last_i.
            last_i = if (last_i >> k) & 1 == 1 {
                last_i - x
            } else {
                last_i + x
            };
            if last_i < n && a[last_i].max > last {
                // Update last accordingly.
                last = a[last_i].max;
            }
            k += 1;
        }
        Some(k - 1)
    }

    /// Lazily yields the indices of all intervals overlapping the query
    /// `[st, en)`, in increasing index (i.e. start) order. Use
    /// [`IITree::start`]/[`IITree::end`]/[`IITree::data`] to look each one up.
    ///
    /// # Panics
    ///
    /// Panics if [`IITree::index`] has not been called since the last
    /// [`IITree::add`].
    pub fn overlap(&self, st: S, en: S) -> Overlap<'_, S, T> {
        assert!(
            self.indexed,
            "IITree::index() must be called before IITree::overlap()"
        );
        let mut it = Overlap {
            tree: self,
            st,
            en,
            stack: [StackCell::default(); 64],
            t: 0,
            i: 0,
            i1: 0,
        };
        if let Some(k) = self.max_level {
            // Push the root; this is a top-down traversal.
            it.push(StackCell {
                x: (1usize << k) - 1,
                k,
                w: false,
            });
        }
        it
    }

    /// `true` if any interval overlaps `[st, en)`.
    pub fn has_overlap(&self, st: S, en: S) -> bool {
        self.overlap(st, en).next().is_some()
    }

    /// Number of intervals in the tree.
    pub fn len(&self) -> usize {
        self.a.len()
    }

    /// `true` if the tree holds no intervals.
    pub fn is_empty(&self) -> bool {
        self.a.is_empty()
    }

    /// Start of interval `i`.
    pub fn start(&self, i: usize) -> S {
        self.a[i].st
    }

    /// End of interval `i`.
    pub fn end(&self, i: usize) -> S {
        self.a[i].en
    }

    /// Payload of interval `i`.
    pub fn data(&self, i: usize) -> &T {
        &self.a[i].data
    }

    /// Mutable payload of interval `i`. Does not invalidate the index (the
    /// coordinates are untouched).
    pub fn data_mut(&mut self, i: usize) -> &mut T {
        &mut self.a[i].data
    }
}

#[derive(Clone, Copy, Default)]
struct StackCell {
    /// Node index (may be `>= a.len()` for a missing node of the complete tree).
    x: usize,
    /// Level of the node.
    k: usize,
    /// `true` once the left child has been processed.
    w: bool,
}

/// Iterator returned by [`IITree::overlap`]. Each call to `next` resumes the
/// top-down traversal from where the previous one left off, so results are
/// produced on demand without an intermediate `Vec`.
pub struct Overlap<'a, S, T> {
    tree: &'a IITree<S, T>,
    st: S,
    en: S,
    stack: [StackCell; 64],
    /// Stack height.
    t: usize,
    /// When `i < i1`, we are inside a small subtree being scanned linearly
    /// over `a[i..i1]`.
    i: usize,
    i1: usize,
}

impl<S: Copy + PartialOrd, T> Overlap<'_, S, T> {
    fn push(&mut self, c: StackCell) {
        self.stack[self.t] = c;
        self.t += 1;
    }
}

impl<S: Copy + PartialOrd, T> Iterator for Overlap<'_, S, T> {
    type Item = usize;

    fn next(&mut self) -> Option<usize> {
        let a = &self.tree.a;
        let n = a.len();
        let (st, en) = (self.st, self.en);
        loop {
            // Finish the linear scan of a small subtree, if one is in progress.
            while self.i < self.i1 && a[self.i].st < en {
                let i = self.i;
                self.i += 1;
                if st < a[i].en {
                    return Some(i);
                }
            }
            self.i1 = self.i; // scan done (or exhausted by the `st < en` cutoff)

            if self.t == 0 {
                return None;
            }
            self.t -= 1;
            let z = self.stack[self.t];
            if z.k <= 3 {
                // We are in a small subtree; traverse every node in it.
                let i0 = z.x >> z.k << z.k;
                let i1 = i0 + (1usize << (z.k + 1)) - 1;
                self.i = i0;
                self.i1 = if i1 >= n { n } else { i1 };
            } else if !z.w {
                // Left child not processed yet.
                let y = z.x - (1usize << (z.k - 1)); // NB: y may be >= n
                // Re-add node z.x, marking the left child as processed.
                self.push(StackCell { w: true, ..z });
                // Push the left child if y is out of range or may overlap the query.
                if y >= n || a[y].max > st {
                    self.push(StackCell {
                        x: y,
                        k: z.k - 1,
                        w: false,
                    });
                }
            } else if z.x < n && a[z.x].st < en {
                // Need to push the right child; it is visited on the next call,
                // after z.x itself has been reported, keeping the output sorted.
                self.push(StackCell {
                    x: z.x + (1usize << (z.k - 1)),
                    k: z.k - 1,
                    w: false,
                });
                if st < a[z.x].en {
                    return Some(z.x);
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Tiny deterministic LCG so the tests need no external crate.
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }
        fn below(&mut self, n: u64) -> u64 {
            self.next() % n
        }
    }

    fn build(rng: &mut Lcg, n: usize, span: u64, max_len: u64) -> IITree<i64, usize> {
        let mut t = IITree::new();
        for i in 0..n {
            let st = rng.below(span) as i64;
            let en = st + 1 + rng.below(max_len) as i64;
            t.add(st, en, i);
        }
        t.index();
        t
    }

    fn naive(t: &IITree<i64, usize>, st: i64, en: i64) -> Vec<usize> {
        (0..t.len())
            .filter(|&i| t.start(i) < en && st < t.end(i))
            .collect()
    }

    #[test]
    fn matches_naive_for_many_sizes() {
        let mut rng = Lcg(42);
        for &n in &[
            0usize, 1, 2, 3, 4, 5, 6, 7, 8, 9, 15, 16, 17, 31, 32, 33, 63, 64, 65, 100, 127, 128,
            129, 255, 256, 257, 1000, 4097,
        ] {
            let t = build(&mut rng, n, 10_000, 500);
            // Indices are in start order after indexing.
            for i in 1..t.len() {
                assert!(t.start(i - 1) <= t.start(i));
            }
            for _ in 0..200 {
                let st = rng.below(10_500) as i64 - 200;
                let en = st + rng.below(800) as i64;
                let got: Vec<usize> = t.overlap(st, en).collect();
                assert_eq!(got, naive(&t, st, en), "n={n} query=[{st},{en})");
                assert_eq!(t.has_overlap(st, en), !got.is_empty());
            }
        }
    }

    #[test]
    fn nested_and_long_intervals() {
        // Long intervals whose max must propagate correctly up the tree.
        let mut rng = Lcg(7);
        for &n in &[10usize, 50, 300, 1023, 1024, 1025] {
            let t = build(&mut rng, n, 1_000_000, 900_000);
            for _ in 0..100 {
                let st = rng.below(1_100_000) as i64;
                let en = st + rng.below(5000) as i64;
                let got: Vec<usize> = t.overlap(st, en).collect();
                assert_eq!(got, naive(&t, st, en), "n={n} query=[{st},{en})");
            }
        }
    }

    #[test]
    fn empty_and_point_queries() {
        let t: IITree<u32, ()> = IITree::new();
        let mut t = t;
        t.index();
        assert!(t.is_empty());
        assert_eq!(t.overlap(0, 100).count(), 0);

        let mut t = IITree::new();
        t.add(10u32, 20, "a");
        t.add(15, 30, "b");
        t.add(40, 50, "c");
        t.index();
        let hits: Vec<&str> = t.overlap(18, 19).map(|i| *t.data(i)).collect();
        assert_eq!(hits, ["a", "b"]);
        // Half-open: a query ending at a start / starting at an end misses it.
        assert_eq!(t.overlap(30, 30).count(), 0);
        assert_eq!(t.overlap(5, 10).count(), 0);
        assert_eq!(t.overlap(20, 20).map(|i| *t.data(i)).collect::<Vec<_>>(), ["b"]);
        assert_eq!(t.overlap(19, 20).map(|i| *t.data(i)).collect::<Vec<_>>(), ["a", "b"]);
        assert_eq!(t.overlap(20, 21).map(|i| *t.data(i)).collect::<Vec<_>>(), ["b"]);
        assert_eq!(t.overlap(30, 40).count(), 0);
        assert_eq!(t.overlap(0, 100).count(), 3);
        // The iterator is lazy: taking one hit is fine.
        assert_eq!(t.overlap(0, 100).next(), Some(0));
        *t.data_mut(2) = "z";
        assert_eq!(*t.data(2), "z");
    }

    #[test]
    #[should_panic(expected = "index()")]
    fn overlap_before_index_panics() {
        let mut t = IITree::new();
        t.add(1i32, 2, ());
        let _ = t.overlap(0, 5);
    }
}
