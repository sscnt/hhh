//
//  GPUFattalFilter.m
//  HDRtest
//
//  Created by SSC on 2014/01/18.
//  Copyright (c) 2014年 SSC. All rights reserved.
//

#import "GPUFattalFilter.h"

@implementation GPUFattalFilter

NSString *const kFattalFragmentShaderString = SHADER_STRING
(
 precision highp float;
 varying vec2 textureCoordinate;
 varying vec2 leftTextureCoordinate;
 varying vec2 rightTextureCoordinate;
 
 varying vec2 topTextureCoordinate;
 varying vec2 topLeftTextureCoordinate;
 varying vec2 topRightTextureCoordinate;
 
 varying vec2 bottomTextureCoordinate;
 varying vec2 bottomLeftTextureCoordinate;
 varying vec2 bottomRightTextureCoordinate;
 
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
     mediump vec3 yuv = rgb2yuv(pixel.rgb);
     
     
     mediump vec3 topPixel = texture2D(inputImageTexture, topTextureCoordinate).rgb;
     mediump vec3 rightPixel = texture2D(inputImageTexture, rightTextureCoordinate).rgb;
     
     mediump float topLuminance = luminance(topPixel);
     mediump float rightLuminance = luminance(rightPixel);
     mediump float luminance = luminance(pixel.rgb);
     
     mediump float topLogLuminance = logLuminance(topPixel);
     mediump float rightLogLuminance = logLuminance(rightPixel);
     mediump float logLuminance = logLuminance(pixel.rgb);
     
     mediump vec2 G;
     G.x = topLogLuminance - logLuminance;
     G.y = rightLogLuminance - logLuminance;
     
     mediump float div = G.x + G.y;
     mediump float absDiv = abs(div);
     mediump float phi = 0.0;
     mediump float alpha = 0.6;
     mediump float beta = 0.8;
     
     if(absDiv == 0.0){
         
     }else{
         phi = (alpha / absDiv) * pow(absDiv / alpha, beta);
     }
     
     div *= phi;
     
     div = (topLogLuminance + rightLogLuminance - div) / 2.0;
     
     div = pow(2.71828182846, div);
     
     
     yuv.x = div;
     

     pixel.rgb = yuv2rgb(yuv);
     
     // Save the result
     gl_FragColor = pixel;
 }
 );


- (id)init;
{
    if (!(self = [super initWithFragmentShaderFromString:kFattalFragmentShaderString]))
    {
        return nil;
    }

    return self;
}


@end
