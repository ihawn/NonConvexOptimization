// <copyright file="BinaryMinHeap.cs">
//  Copyright (c) Eric Force.
//  Licensed under the MIT license. See LICENSE file in the project root for full license information.
// </copyright>
// <summary>
//  This class implements the binary heap algorithm where the min value is always the root.
// </summary>

namespace BinaryHeap
{
    using System;
    using System.Collections.Generic;
    using System.Diagnostics.CodeAnalysis;
    using System.Runtime.InteropServices;
    using System.Text;

    /// <summary>
    /// Implements the binary heap algorithm where the min value is always the root. A binary heap allows on average O(1) insertion while representing the tree as a single dimensional List.
    /// </summary>
    /// <typeparam name="T">The type of elements in the heap.</typeparam>
    [ComVisible(true)]
    [Guid("05AC1509-63C5-4922-8B3E-364E6797F1C1")]
    [ClassInterface(ClassInterfaceType.None)]
    public class BinaryMinHeap<T> : IBinaryHeap<T>
        where T : IComparable<T>
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="BinaryMinHeap{T}"/> class with an empty heap.
        /// </summary>
        public BinaryMinHeap()
        {
            this.Heap = new List<T>();
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="BinaryMinHeap{T}"/> class that contains elements copied from the specified collection.
        /// <para>The heap is built in worse case O(nlog(n)) time.</para>
        /// </summary>
        /// <param name="collection">The collection whose elements are copied to the new list.</param>
        /// <exception cref="ArgumentNullException">collection is null.</exception>
        public BinaryMinHeap(IEnumerable<T> collection)
        {
            if (collection == null)
            {
                throw new ArgumentNullException(nameof(collection));
            }

            this.Heap = new List<T>();

            foreach (T item in collection)
            {
                this.Insert(item);
            }
        }

        /// <summary>
        /// Initializes a new instance of the <see cref="BinaryMinHeap{T}"/> class that contains elements copied from the specified collection.
        /// <para>The heap is built in worse case O(n) time (Floyd).</para>
        /// </summary>
        /// <param name="collection">The collection whose elements are copied to the new list.</param>
        /// <exception cref="ArgumentNullException">collection is null.</exception>
        public BinaryMinHeap(ICollection<T> collection)
        {
            if (collection == null)
            {
                throw new ArgumentNullException(nameof(collection));
            }

            this.Heap = new List<T>(collection);

            for (int i = FindParentIndex(this.Heap.Count); i >= 0; i--)
            {
                this.DownHeap(i);
            }
        }

        /// <summary>
        /// Gets the number of elements contained in the <see cref="BinaryMinHeap{T}"/>.
        /// </summary>
        public int Count => this.Heap.Count;

        /// <summary>
        /// Gets the inner <see cref="IList{T}"/> for the heap. The binary tree is represented as a single dimensional array.
        /// </summary>
        protected IList<T> Heap { get; private set; }

        /// <summary>
        /// Inserts an item into the <see cref="BinaryMinHeap{T}"/>.
        /// </summary>
        /// <param name="item">The object to add to the <see cref="BinaryMinHeap{T}"/>.</param>
        public void Insert(T item)
        {
            this.Heap.Add(item);

            // move the last inserted item up the heap if needed
            this.UpHeap(this.Heap.Count - 1);
        }

        /// <summary>
        /// Returns and removes the root or lowest value in the heap.
        /// </summary>
        /// <returns>The root or lowest value in the heap.</returns>
        /// <exception cref="InvalidOperationException">The <see cref="BinaryMinHeap{T}"/> is empty.</exception>
        public T Extract()
        {
            if (this.Count == 0)
            {
                throw new InvalidOperationException("Empty heap.");
            }

            // swap the root node with the last node and remove the root before returning it
            T root = this.Heap[0];
            this.Heap[0] = this.Heap[this.Heap.Count - 1];
            this.Heap.RemoveAt(this.Heap.Count - 1);

            if (this.Heap.Count == 0)
            {
                // no more elements in heap, return root and don't down heap
                return root;
            }

            // move the new root down if needed
            this.DownHeap(0);

            return root;
        }

        /// <summary>
        /// Returns the object at the root of the <see cref="BinaryMinHeap{T}"/> without removing it.
        /// </summary>
        /// <returns>The object at the root of the <see cref="BinaryMinHeap{T}"/>.</returns>
        public T Peek()
        {
            if (this.Count == 0)
            {
                throw new InvalidOperationException("Empty heap.");
            }

            return this.Heap[0];
        }

        /// <summary>
        /// Returns a string that represents the current <see cref="BinaryMinHeap{T}"/> instance.
        /// </summary>
        /// <returns>A string that represents the current <see cref="BinaryMinHeap{T}"/> instance.</returns>
        public override string ToString()
        {
            StringBuilder sb = new StringBuilder();

            for (int i = 0; i < this.Heap.Count; ++i)
            {
                sb.Append(this.Heap[i]).Append(", ");
            }

            sb.Append("count = ").Append(this.Heap.Count);

            return sb.ToString();
        }

        /// <summary>
        /// Iteratively swaps the current node and the child node if the parent node isn't the smallest between the two children.
        /// </summary>
        /// <param name="currentIndex">The starting node's index.</param>
        protected void UpHeap(int currentIndex)
        {
            if (currentIndex == 0)
            {
                // at root
                return;
            }

            // find the parent of the node and its child
            int parentIndex = FindParentIndex(currentIndex);

            T current = this.Heap[currentIndex];
            T parent = this.Heap[parentIndex];

            if (current.CompareTo(parent) < 0)
            {
                // child is less than parent so swap the parent and child and move the new parent up if needed
                this.Heap[parentIndex] = current;
                this.Heap[currentIndex] = parent;

                this.UpHeap(parentIndex);
            }
        }

        /// <summary>
        /// Iteratively swaps the current node and the largest child node if the current node isn't the smallest between its two children.
        /// </summary>
        /// <param name="currentIndex">The starting node's index.</param>
        protected void DownHeap(int currentIndex)
        {
            int rightIndex = FindRightChildIndex(currentIndex);
            int leftIndex = FindLeftChildIndex(currentIndex);

            T current = this.Heap[currentIndex];
            T right = rightIndex < this.Heap.Count ? this.Heap[rightIndex] : default(T); // check if a right node exists
            T left = leftIndex < this.Heap.Count ? this.Heap[leftIndex] : default(T); // check if a left node exists

            // only do a compare if a node actually exists. If a node doesn't exist, it's effectively lower than the sibling node or current working node
            if ((rightIndex >= this.Heap.Count || LessThanOrEqual(current, right)) && (leftIndex >= this.Heap.Count || LessThanOrEqual(current, left)))
            {
                // the parent is the smallest of the three so the heap is valid
                return;
            }
            else if ((rightIndex < this.Heap.Count && LessThanOrEqual(right, current)) && LessThanOrEqual(right, left))
            {
                // the left is assumed to exist if the first if is false since the first check in this if is for the right node
                // the right is smaller so swap the right with the parent and move the new right down if needed
                this.Heap[currentIndex] = right;
                this.Heap[rightIndex] = current;

                this.DownHeap(rightIndex);
            }
            else
            {
                // the left is smaller so swap the left with the parent and move the new left down if needed
                this.Heap[currentIndex] = left;
                this.Heap[leftIndex] = current;

                this.DownHeap(leftIndex);
            }
        }

        // static helpers
        protected static bool LessThanOrEqual(T first, T second)
        {
            return first.CompareTo(second) <= 0;
        }

        protected static int FindLeftChildIndex(int index)
        {
            return (2 * index) + 1;
        }

        protected static int FindRightChildIndex(int index)
        {
            return (2 * index) + 2;
        }

        protected static int FindParentIndex(int index)
        {
            return (index - 1) / 2;
        }
    }
}
