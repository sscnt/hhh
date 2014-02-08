//
//  GPUDefaultFilter.h
//  HDRtest
//
//  Created by SSC on 2014/01/19.
//  Copyright (c) 2014年 SSC. All rights reserved.
//

#import "GPUImageTwoInputCrossTextureSamplingFilter.h"
#import "GPUImageFilterGroup.h"

extern NSString *const kGPUGaussSeidelFilterFragmentShaderString;

@interface GPUGaussSeidelFilter : GPUImageTwoInputCrossTextureSamplingFilter
{
    GLint mixUniform;
    
    GLuint secondFilterOutputTexture;
    GLuint secondFilterFramebuffer;
}

// The number of times to propagate the gradients.
// Crank this up to 100 or even 1000 if you want to get anywhere near convergence.  Yes, this will be slow.
@property(readwrite, nonatomic) NSUInteger numIterations;

@end
