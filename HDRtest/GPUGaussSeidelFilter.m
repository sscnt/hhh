//
//  GPUDefaultFilter.m
//  HDRtest
//
//  Created by SSC on 2014/01/19.
//  Copyright (c) 2014年 SSC. All rights reserved.
//

#import "GPUGaussSeidelFilter.h"

@implementation GPUGaussSeidelFilter

NSString *const kDefaultFragmentShaderString = SHADER_STRING
(
 precision highp float;
 varying vec2 textureCoordinate;
 uniform sampler2D inputImageTexture;
 uniform mediump float mtplML;
 
 mediump float luminance(mediump vec3 pixel)
{
    return 0.299 * pixel.r + 0.587 * pixel.g + 0.114 * pixel.b;
}
 
 mediump float logLuminance(mediump vec3 pixel){
     return log(0.00001 + luminance(pixel));
 }
 
 
 vec3 rgb2yuv(vec3 rgb){
     mediump float r = rgb.r * 0.8588 + 0.0625;
     mediump float g = rgb.g * 0.8588 + 0.0625;
     mediump float b = rgb.b * 0.8588 + 0.0625;
     
     mediump float y = 0.299 * r + 0.587 * g + 0.114 * b;
     mediump float u = -0.169 * r - 0.331 * g + 0.500 * b;
     mediump float v = 0.500 * r - 0.419 * g - 0.081 * b;
     return vec3(y, u, v);
 }
 
 vec3 yuv2rgb(vec3 yuv){
     mediump float r = yuv.x + 1.402 * yuv.z - 0.0625;
     mediump float g = yuv.x - 0.344 * yuv.y - 0.714 * yuv.z - 0.0625;
     mediump float b = yuv.x + 1.772 * yuv.y - 0.0625;
     
     r *= 1.164;
     g *= 1.164;
     b *= 1.164;
     
     r = max(0.0, min(1.0, r));
     g = max(0.0, min(1.0, g));
     b = max(0.0, min(1.0, b));
     return vec3(r, g, b);
 }
 
 void main()
 {
     mediump vec4 pixel = texture2D(inputImageTexture, textureCoordinate);
     mediump vec3 yuv;
     mediump float alpha;
     
     mediump vec3 multiplied = pixel.rgb * pixel.rgb;
     yuv = rgb2yuv(multiplied);
     alpha = 1.0 / (1.0 + exp(-(yuv.x * 12.0 - 6.0)));
     alpha = max(min(1.0, alpha), 0.0);
     
     multiplied = multiplied * alpha + pixel.rgb * (1.0 - alpha);
     
     mediump vec3 screened = 1.0 - ((1.0 - pixel.rgb) * (1.0 - pixel.rgb));
     yuv = rgb2yuv(screened);
     alpha = 1.0 / (1.0 + exp((yuv.x * 12.0 - 6.0)));
     alpha = max(min(1.0, alpha), 0.0);
     
     screened = screened * alpha + pixel.rgb * (1.0 - alpha);
     pixel.rgb = screened;
     
     yuv.x = 0.5;
     pixel.rgb = yuv2rgb(yuv);
     
     // Save the result
     gl_FragColor = pixel;
 }
 );


- (id)init;
{
    if (!(self = [super initWithFragmentShaderFromString:kDefaultFragmentShaderString]))
    {
        return nil;
    }
    
    return self;
}
@end
