// <copyright file="IBinaryHeap.cs">
//  Copyright (c) Eric Force.
//  Licensed under the MIT license. See LICENSE file in the project root for full license information.
// </copyright>
// <summary>
//  This interface defines the members used for a binary heap algorithm.
// </summary>

namespace BinaryHeap
{
    using System;
    using System.Runtime.InteropServices;

    /// <summary>
    /// Represents a binary heap algorithm. A binary heap allows on average O(1) insertion while representing the tree as a single dimensional List.
    /// </summary>
    /// <typeparam name="T">The type of elements in the heap.</typeparam>
    [ComVisible(true)]
    [Guid("C786D982-B693-45F0-953D-56BD98D5D1BE")]
    [InterfaceType(ComInterfaceType.InterfaceIsDual)]
    public interface IBinaryHeap<T>
    {
        /// <summary>
        /// Gets the number of elements contained in the <see cref="IBinaryHeap{T}"/>.
        /// </summary>
        int Count { get; }

        /// <summary>
        /// Inserts an item into the <see cref="IBinaryHeap{T}"/>.
        /// </summary>
        /// <param name="item">The object to add to the <see cref="IBinaryHeap{T}"/>.</param>
        void Insert(T item);

        /// <summary>
        /// Returns and removes the root or the heap.
        /// </summary>
        /// <returns>The root in the heap.</returns>
        /// <exception cref="InvalidOperationException">The <see cref="IBinaryHeap{T}"/> is empty.</exception>
        T Extract();

        /// <summary>
        /// Returns the object at the root of the <see cref="IBinaryHeap{T}"/> without removing it.
        /// </summary>
        /// <returns>The object at the root of the <see cref="IBinaryHeap{T}"/>.</returns>
        T Peek();
    }
}
