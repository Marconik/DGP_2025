/***************************************************************************
 # Copyright (c) 2015-24, NVIDIA CORPORATION. All rights reserved.
 #
 # Redistribution and use in source and binary forms, with or without
 # modification, are permitted provided that the following conditions
 # are met:
 #  * Redistributions of source code must retain the above copyright
 #    notice, this list of conditions and the following disclaimer.
 #  * Redistributions in binary form must reproduce the above copyright
 #    notice, this list of conditions and the following disclaimer in the
 #    documentation and/or other materials provided with the distribution.
 #  * Neither the name of NVIDIA CORPORATION nor the names of its
 #    contributors may be used to endorse or promote products derived
 #    from this software without specific prior written permission.
 #
 # THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS "AS IS" AND ANY
 # EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 # IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 # PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 # CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 # EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 # PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 # PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 # OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 # (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 # OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 **************************************************************************/
#pragma once
#include <functional>
#include <memory>
#include <vector>

#include "Core/Macros.h"
#include "LightBVHTypes.slang"
#include "Scene/Lights/LightCollection.h"
#include "utils/Math/AABB.h"
#include "utils/Math/Vector.h"

namespace USTC_CG {
class LightBVHBuilder;

/** Utility class representing a light sampling BVH.

    This is binary BVH over all emissive triangles as described by Moreau and
   Clarberg, "Importance Sampling of Many Lights on the GPU", Ray Tracing Gems,
   Ch. 18, 2019.

    Before being used, the BVH needs to have been built using
   LightBVHBuilder::build(). The data can be both used on the CPU (using
   traverseBVH() or on the GPU by:
      1. import LightBVH;
      2. Declare a variable of type LightBVH in your shader.
      3. Call bindShaderData() to bind the BVH resources.

    TODO: Rename all things 'triangle' to 'light' as the BVH can be used for
   other light types.
*/
class HD_USTC_CG_API LightBVH {
   public:
    struct NodeLocation {
        uint32_t nodeIndex;
        uint32_t depth;

        NodeLocation() : nodeIndex(0), depth(0)
        {
        }
        NodeLocation(uint32_t _nodeIndex, uint32_t _depth)
            : nodeIndex(_nodeIndex),
              depth(_depth)
        {
        }
    };

    /** Function called on each node by traverseBVH().
        \param[in] location The location of the node in the tree.
        \return True if the traversal should continue, false otherwise.
    */
    using NodeFunction = std::function<bool(const NodeLocation& location)>;

    /** Constructor.
        \param[in] pDevice GPU device.
        \param[in] pLightCollection The light collection around which the BVH
       will be built.
    */
    LightBVH(
        nvrhi::DeviceHandle pDevice,
        const ref<const ILightCollection>& pLightCollection);

    /** Returns the LightCollection the LightBVH is built on.
     */
    ref<const ILightCollection> getLightCollection() const
    {
        return mpLightCollection;
    }

    /** Refit all the BVH nodes to the underlying geometry, without changing the
       hierarchy. The BVH needs to have been built before trying to refit it.
        \param[in] pRenderContext The render context.
    */
    void refit(RenderContext* pRenderContext);

    /** Perform a depth-first traversal of the BVH and run a function on each
       node. \param[in] evalInternal Function called on each internal node.
        \param[in] evalLeaf Function called on each leaf node.
        \param[in] rootNodeIndex The index of the node to start traversing.
    */
    void traverseBVH(
        const NodeFunction& evalInternal,
        const NodeFunction& evalLeaf,
        uint32_t rootNodeIndex = 0);

    struct BVHStats {
        std::vector<uint32_t>
            nodeCountPerLevel;  ///< For each level in the tree, how many nodes
                                ///< are there.
        std::vector<uint32_t>
            leafCountPerTriangleCount;  ///< For each amount of triangles, how
                                        ///< many leaf nodes contain that many
                                        ///< triangles.

        uint32_t treeHeight = 0;  ///< Number of edges on the longest path
                                  ///< between the root node and a leaf.
        uint32_t minDepth = 0;    ///< Number of edges on the shortest path
                                  ///< between the root node and a leaf.
        uint32_t byteSize = 0;    ///< Number of bytes occupied by the BVH.
        uint32_t internalNodeCount =
            0;  ///< Number of internal nodes inside the BVH.
        uint32_t leafNodeCount = 0;  ///< Number of leaf nodes inside the BVH.
        uint32_t triangleCount = 0;  ///< Number of triangles inside the BVH.
    };

    /** Returns stats.
     */
    const BVHStats& getStats() const
    {
        return mBVHStats;
    }

    /** Is the BVH valid.
        \return true if the BVH is ready for use.
    */
    bool isValid() const
    {
        return mIsValid;
    }

    /** Bind the light BVH into a shader variable.
        \param[in] var The shader variable to set the data into.
    */
    void bindShaderData(const ShaderVar& var) const;

   protected:
    void finalize();
    void computeStats();
    void updateNodeIndices();

    void uploadCPUBuffers(
        const std::vector<uint32_t>& triangleIndices,
        const std::vector<uint64_t>& triangleBitmasks);
    void syncDataToCPU() const;

    /** Invalidate the BVH.
     */
    void clear();

    struct RefitEntryInfo {
        uint32_t offset = 0;  ///< Offset into the 'mpNodeIndicesBuffer' buffer.
        uint32_t count = 0;   ///< The number of nodes at each level.
    };

    // Internal state
    nvrhi::DeviceHandle mpDevice;
    ref<const ILightCollection> mpLightCollection;

    ref<ComputePass>
        mLeafUpdater;  ///< Compute pass for refitting the leaf nodes.
    ref<ComputePass>
        mInternalUpdater;  ///< Compute pass for refitting internal nodes.

    // CPU resources
    mutable std::vector<PackedNode>
        mNodes;  ///< CPU-side copy of packed BVH nodes.
    std::vector<uint32_t>
        mNodeIndices;  ///< Array of all node indices sorted by tree depth.
    std::vector<RefitEntryInfo>
        mPerDepthRefitEntryInfo;  ///< Array containing for each level the
                                  ///< number of internal nodes as well as the
                                  ///< corresponding offset into
                                  ///< 'mpNodeIndicesBuffer'; the very last
                                  ///< entry contains the same data, but for all
                                  ///< leaf nodes instead.
    uint32_t mMaxTriangleCountPerLeaf =
        0;  ///< After the BVH is built, this contains the maximum light count
            ///< per leaf node.
    BVHStats mBVHStats;
    bool mIsValid = false;  ///< True when the BVH has been built.
    mutable bool mIsCpuDataValid = false;  ///< Indicates whether the CPU-side
                                           ///< data matches the GPU buffers.

    // GPU resources
    nvrhi::BufferHandle mpBVHNodesBuffer;  ///< Buffer holding all BVH nodes.
    nvrhi::BufferHandle
        mpTriangleIndicesBuffer;  ///< Triangle indices sorted by leaf node.
                                  ///< Each leaf node refers to a contiguous
                                  ///< array of triangle indices.
    nvrhi::BufferHandle
        mpTriangleBitmasksBuffer;  ///< Array containing the per triangle bit
                                   ///< pattern retracing the tree traversal to
                                   ///< reach the triangle: 0=left child,
                                   ///< 1=right child.
    nvrhi::BufferHandle
        mpNodeIndicesBuffer;  ///< Buffer holding all node indices sorted by
                              ///< tree depth. This is used for BVH refit.

    friend LightBVHBuilder;
};
}  // namespace USTC_CG
