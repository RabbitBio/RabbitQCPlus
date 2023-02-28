#ifndef CARE_CPU_BEST_ALIGNMENT_HPP
#define CARE_CPU_BEST_ALIGNMENT_HPP

#include <config.hpp>
#include <hpc_helpers.cuh>

namespace care{


    enum class AlignmentOrientation : char {Forward=1, ReverseComplement=2, None=3};

    HOSTDEVICEQUALIFIER
    __inline__
    AlignmentOrientation chooseBestAlignmentOrientation(int fwd_alignment_overlap,
    			int revc_alignment_overlap,
    			int fwd_alignment_nops,
    			int revc_alignment_nops,
    			bool fwd_alignment_isvalid,
    			bool revc_alignment_isvalid,
    			int anchorLength,
    			int /*querylength*/,
    			float min_overlap_ratio,
    			int min_overlap,
    			float maxErrorRate){

    	AlignmentOrientation retval = AlignmentOrientation::None;

    	const int minimumOverlap = int(anchorLength * min_overlap_ratio) > min_overlap
    	                           ? int(anchorLength * min_overlap_ratio) : min_overlap;

    	// choose alignment with smallest error rate in overlap and overlaplength >= minimumOverlap and error rate in overlap < maxErrorRate

    	if(fwd_alignment_isvalid && fwd_alignment_overlap >= minimumOverlap) {
    		if(revc_alignment_isvalid && revc_alignment_overlap >= minimumOverlap) {
    			const float ratio = (float)fwd_alignment_nops / fwd_alignment_overlap;
    			const float revcomplratio = (float)revc_alignment_nops / revc_alignment_overlap;

    			if(ratio < revcomplratio) {
    				if(ratio < maxErrorRate) {
    					retval = AlignmentOrientation::Forward;
    				}
    			}else if(revcomplratio < ratio) {
    				if(revcomplratio < maxErrorRate) {
    					retval = AlignmentOrientation::ReverseComplement;
    				}
    			}else{
    				if(ratio < maxErrorRate) {
    					// both have same mismatch ratio, choose longest overlap
    					if(fwd_alignment_overlap > revc_alignment_overlap) {
    						retval = AlignmentOrientation::Forward;
    					}else{
    						retval = AlignmentOrientation::ReverseComplement;
    					}
    				}
    			}
    		}else{
    			if((float)fwd_alignment_nops / fwd_alignment_overlap < maxErrorRate) {
    				retval = AlignmentOrientation::Forward;
    			}
    		}
    	}else{
    		if(revc_alignment_isvalid && revc_alignment_overlap >= minimumOverlap) {
    			if((float)revc_alignment_nops / revc_alignment_overlap < maxErrorRate) {
    				retval = AlignmentOrientation::ReverseComplement;
    			}
    		}
    	}

    	return retval;
    }

    template<class Alignment>
    AlignmentOrientation chooseBestAlignmentOrientation(const Alignment& fwdAlignment,
                                        const Alignment& revcmplAlignment,
                                        int anchorLength,
                                        int querylength,
                                        float min_overlap_ratio,
                                        int min_overlap,
                                        float estimatedErrorrate){

        return chooseBestAlignmentOrientation(fwdAlignment.get_overlap(),
        			revcmplAlignment.get_overlap(),
        			fwdAlignment.get_nOps(),
        			revcmplAlignment.get_nOps(),
        			fwdAlignment.get_isValid(),
        			revcmplAlignment.get_isValid(),
        			anchorLength,
        			querylength,
        			min_overlap_ratio,
        			min_overlap,
        			estimatedErrorrate * 4.0f);
    }

}



#endif
