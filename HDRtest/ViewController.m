//
//  ViewController.m
//  HDRtest
//
//  Created by SSC on 2014/01/18.
//  Copyright (c) 2014年 SSC. All rights reserved.
//

#import "ViewController.h"'
#import "GPUImageGaussianBlurFilter.h"

@interface ViewController ()

@end

@implementation ViewController

- (void)viewDidLoad
{
    [super viewDidLoad];
	// Do any additional setup after loading the view, typically from a nib.
    
    UIButton* button = [[UIButton alloc] initWithFrame:CGRectMake(0.0f, 0.0f, 320.0f, 44.0f)];
    [button setTitle:@"Open" forState:UIControlStateNormal];
    [button setTitleColor:[UIColor blackColor] forState:UIControlStateNormal];
    [button addTarget:self action:@selector(didButtonPress) forControlEvents:UIControlEventTouchUpInside];
    [self.view addSubview:button];
}

- (void)didButtonPress
{
    UIImagePickerController *pickerController = [[UIImagePickerController alloc] init];
    pickerController.delegate = self;
    [pickerController setSourceType:UIImagePickerControllerSourceTypePhotoLibrary];
    pickerController.modalTransitionStyle = UIModalTransitionStyleCrossDissolve;
    [self presentViewController:pickerController animated:YES completion:nil];
}

#pragma mark delegates

- (void)imagePickerController:(UIImagePickerController *)picker didFinishPickingMediaWithInfo:(NSDictionary *)info
{
    [picker dismissViewControllerAnimated:YES completion:nil];
    self.loadedImage = [info objectForKey:UIImagePickerControllerOriginalImage];
    if(!self.loadedImage){
        NSURL* imageurl = [info objectForKey:UIImagePickerControllerReferenceURL];
        ALAssetsLibrary* library = [[ALAssetsLibrary alloc] init];
        __block ViewController* _self = self;
        [library assetForURL:imageurl
                 resultBlock: ^(ALAsset *asset)
         {
             ALAssetRepresentation *representation;
             representation = [asset defaultRepresentation];
             _self.loadedImage = [[UIImage alloc] initWithCGImage:representation.fullResolutionImage];
             [_self didLoadImage];
         }
                failureBlock:^(NSError *error)
         {
             NSLog(@"error:%@", error);
         }];
    } else {
        if(picker.sourceType == UIImagePickerControllerSourceTypeCamera){
            UIImageWriteToSavedPhotosAlbum(self.loadedImage, nil, nil, nil);
        }
        //[picker.delegate performSelector:@selector(prepareForEditor) withObject:nil];
        [self didLoadImage];
    }
}

- (void)didLoadImage
{
    UIImage* result;
    
    if(true){
        result = [self hdrWithInputImage:self.loadedImage];

        
    }    else{
        

    }
    
    CGFloat base_width = 320.0f;
    
    UIImageView* view = [[UIImageView alloc] initWithImage:self.loadedImage];
    CGFloat height = self.loadedImage.size.height * base_width / self.loadedImage.size.width;
    [view setFrame:CGRectMake(0.0f, 44.0f, base_width, height)];
    [self.view addSubview:view];
    
    NSLog(@"%f, %f", height, self.loadedImage.size.width);
    
    CGFloat y = height + 10.0f;
    
    view = [[UIImageView alloc] initWithImage:result];
    height = self.loadedImage.size.height * base_width / self.loadedImage.size.width;
    [view setFrame:CGRectMake(0.0f, 44.0f + y, base_width, height)];
    [self.view addSubview:view];
}

float rgb2luminance(float red, float green, float blue)
{
    return 0.299 * red + 0.587 * green + 0.114 * blue;
}

float pixel2luminance(unsigned char* pixel)
{
    return rgb2luminance(*(pixel), *(pixel + 1), *(pixel + 2));
}



void multiply(float* r, float* g, float* b)
{
    *r = *r * *r;
    *g = *g * *g;
    *b = *b * *b;
}

void screen2(float* r, float* g, float* b)
{

    *r = 1.0 - (1.0 - *r) * (1.0 - *r);
    *g = 1.0 - (1.0 - *g) * (1.0 - *g);
    *b = 1.0 - (1.0 - *b) * (1.0 - *b);
}

float _log(float v)
{
    return log(v);
    float t = v - 1.0;
    return t - t * t * 0.5 + t * t * t * 0.33333333333333 - t * t * t * t * 0.25;
}

void screen(float* value)
{
    
    float weight = _log(*value);
    weight = weight * weight * 0.05;
    weight = MIN(weight, 1.0);
    weight = 1.0 - weight;
    *value = 1.0 - (1.0 - *value * weight) * (1.0 - *value * weight);
}

void screen_rgb(float* r, float* g, float* b)
{
    float weight = _log(*r);
    weight = weight * weight * 0.05;
    weight = MIN(weight, 1.0);
    weight = 1.0 - weight;
    *r = 1.0 - (1.0 - *r * weight) * (1.0 - *r * weight);
    weight = _log(*g);
    weight = weight * weight * 0.05;
    weight = MIN(weight, 1.0);
    weight = 1.0 - weight;
    *g = 1.0 - (1.0 - *g * weight) * (1.0 - *g * weight);
    weight = _log(*b);
    weight = weight * weight * 0.05;
    weight = MIN(weight, 1.0);
    weight = 1.0 - weight;
    *b = 1.0 - (1.0 - *b * weight) * (1.0 - *b * weight);
}


void unko(unsigned char* rawpixel)
{
    float r, g, b;
    r = (float)*(rawpixel + 0) * 0.00392156862745;
    g = (float)*(rawpixel + 1) * 0.00392156862745;
    b = (float)*(rawpixel + 2) * 0.00392156862745;
    
    screen_rgb(&r, &g, &b);
    *(rawpixel + 0) = MAX(MIN((int)(r * 255.0), 255), 0);
    *(rawpixel + 1) = MAX(MIN((int)(g * 255.0), 255), 0);
    *(rawpixel + 2) = MAX(MIN((int)(b * 255.0), 255), 0);
}

typedef struct
{
    float r;
    float g;
    float b;
} RGB;

typedef struct
{
    float r;
    float g;
    float b;
    float a;
} RGBA;

typedef struct
{
    float y;
    float u;
    float v;
} YUV;

typedef struct
{
    int x;
    int y;
    float ev0;
    float ev2;
    float ev4;
    float ev_2;
    float ev_4;
} Position;

typedef struct
{
    float luminance;
    float luminance_multiplied;
    float luminance_screened;
    float delta_lum_multiplied;
    float delta_lum_screened;
} Pixel;

YUV rgb2yuv(float r, float g, float b)
{
    YUV color = {0.0, 0.0, 0.0};
    color.y = 0.299 * r + 0.587 * g + 0.114 * b;
    color.u = -0.169 * r - 0.331 * g + 0.500 * b;
    color.v = 0.500 * r - 0.419 * g - 0.081 * b;
    return color;
}

RGB yuv2rgb(float y, float u, float v)
{
    RGB color = {0.0, 0.0, 0.0};
    color.r = y + 1.402 * v;
    color.g = y - 0.344 * u - 0.714 * v;
    color.b = y + 1.772 * u;
    return color;
}

void calc(unsigned char* rawpixel, Pixel* pixel)
{
    float r, g, b, _r, _g, _b;
    r = (float)*(rawpixel + 0) * 0.00392156862745;
    g = (float)*(rawpixel + 1) * 0.00392156862745;
    b = (float)*(rawpixel + 2) * 0.00392156862745;
    (*pixel).luminance = rgb2luminance(r, g, b);
    
    _r = r;
    _g = g;
    _b = b;
    multiply(&_r, &_g, &_b);
    (*pixel).luminance_multiplied = rgb2luminance(_r, _g, _b);
    (*pixel).delta_lum_multiplied = (*pixel).luminance - (*pixel).luminance_multiplied;
    
    _r = r;
    _g = g;
    _b = b;
    screen_rgb(&_r, &_g, &_b);
    (*pixel).luminance_screened = rgb2luminance(_r, _g, _b);
    (*pixel).delta_lum_screened = (*pixel).luminance_screened - (*pixel).luminance;
}

void inverse(float input[][100])
{
    float inv_a[100][100]; //ここに逆行列が入る
    float buf; //一時的なデータを蓄える
    int i,j,k; //カウンタ
    int n=100;  //配列の次数
    
    //単位行列を作る
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            inv_a[i][j]=(i==j)?1.0:0.0;
        }
    }
    //掃き出し法
    for(i=0;i<n;i++){
        buf=1/input[i][i];
        for(j=0;j<n;j++){
            input[i][j]*=buf;
            inv_a[i][j]*=buf;
        }
        for(j=0;j<n;j++){
            if(i!=j){
                buf=input[j][i];
                for(k=0;k<n;k++){
                    input[j][k]-=input[i][k]*buf;
                    inv_a[j][k]-=inv_a[i][k]*buf;
                }
            }
        }
    }
}

float absd(float value)
{
    if(value  < 0.0){
        return -value;
    }
    return value;
}

void gradient(float* g_div, float* radiances, float* phis, int width, int height){
    
    float phi, lum;
    
    for (int j = 0; j < height; j++) {
        g_div[j * width + 0] = 0.0;
        g_div[j * width + width - 1] = 0.0;
    }
    
    for (int i = 0; i < width; i++) {
        g_div[0 * width + i] = 0.0;
        g_div[(height - 1) * width + width - 1] = 0.0;
    }
    
    float tmp;
    for (int j = 1 ; j < height - 1; j++)
    {
        for (int i = 1; i < width - 1; i++)
        {
        
            lum = radiances[(j * width) * 4 + i * 4 + 3];

            phi = phis[j * width + i];
            phi = 1.0;
            
            g_div[j * width + i] += phi * (radiances[j * width * 4 + (i - 1) * 4 + 3] + radiances[j * width * 4 + (i + 1) * 4 + 3] - 2.0 * lum);
            g_div[j * width + i] += phi * (radiances[(j - 1) * width * 4 + i * 4 + 3] + radiances[(j + 1) * width * 4 + i * 4 + 3] - 2.0 * lum);
            
            tmp = g_div[j * width + i];
        }
    }
}

void gaussb_(float* result, float* g_div, int width, int height){
    
    bool cflag = true;
    bool* flags = (bool*)malloc(sizeof(bool) * width * height);
    int updated = 1;
    float prev_i, max_i, current_err = 0.0, threshold = 0.0001;
    
    while (updated > 0) {
        updated = 0;
        for (int j = 1 ; j < height - 1; j++)
        {
            for (int i = 1; i < width - 1; i++)
            {
                if(flags[j * width + i] == cflag){
                    continue;
                }
                
                prev_i = result[j * width + i];
                result[j * width + i] = 0.25 * (result[j * width + i + 1] + result[j * width + i - 1] + result[(j + 1) * width + i] + result[(j - 1) * width + i] - g_div[j * width + i]);
                //result[j * width + i] = prev_i + 1.25 * (result[j * width + i] - prev_i);
                if(absd(result[j * width + i]) > max_i){
                    max_i = absd(result[j * width + i]);
                }
                //tmp = result[(j + 1) * width + i] + result[(j - 1) * width + i] + result[j * width + (i + 1)] + result[j * width + (i - 1)] - 4.0 * result[j * width + i];
                current_err = absd(prev_i - result[j * width + i]) / max_i;
                if(current_err > threshold){
                    updated++;
                    flags[j * width + i] = !cflag;
                    flags[j * width + i - 1] = !cflag;
                    flags[j * width + i + 1] = !cflag;
                    flags[(j + 1) * width + i] = !cflag;
                    flags[(j - 1) * width + i] = !cflag;
                    flags[(j + 1) * width + i + 1] = !cflag;
                    flags[(j + 1) * width + i - 1] = !cflag;
                    flags[(j - 1) * width + i + 1] = !cflag;
                    flags[(j - 1) * width + i - 1] = !cflag;
                }else{
                    flags[j * width + i] = cflag;
                }
            }
        }
        NSLog(@"Updated:%d", updated);
    }
    
    free(flags);
}

void expall(float* result, float* radiances, int width, int height){
    float tmp;
    for (int j = 1 ; j < height - 1; j++)
    {
        for (int i = 1; i < width - 1; i++)
        {
            result[j * width + i] = exp(result[j * width + i]);
            radiances[j * width * 4 + i * 4 + 3] = exp(radiances[j * width * 4 + i * 4 + 3]);
            tmp = result[j * width + i] / (1.0 + result[j * width + i]);
            result[j * width + i] = tmp;

        }
    }
}

void calcphis(float* phis, float* radiances, int width, int height, float alpha, float beta){
    float phi, h_x, h_y, tmp, h;
    for (int j = 1 ; j < height - 1; j++)
    {
        for (int i = 1; i < width - 1; i++)
        {
            phi = 1.0;
            h_x = (radiances[j * width * 4 + (i + 1) * 4 + 3] - radiances[j * width * 4 + (i - 1) * 4 + 3]) / 2.0;
            h_y = (radiances[(j + 1) * width * 4 + i * 4 + 3] - radiances[(j - 1) * width * 4 + i * 4 + 3]) / 2.0;
            h = sqrt(h_x * h_x + h_y * h_y);
            
            tmp = (alpha / h) * pow(h / alpha, beta);
            if(tmp > 1.0){
                tmp = 1.0;
            }
            phi *= tmp;
            if(isnan(phi)){
                phi = 1.0;
            }
            phi = MAX(phi, 1.0);
            phis[j * width + i] = phi;
        }
    }

}

RGBA float2RGBA(float input){
    float value = absd(input) / 2.0;
    float tmp;
    RGBA rgba = {0, 0, 0, 0};
    //value *= max;
    //value = floor(value);
    tmp = floor(value / 0.00390625);
    
    rgba.a = tmp / 255.0;
    if(input < 0.0){
        rgba.a += 128.0/255.0;
    }
    value -= tmp * 0.00390625;
    
    tmp = floor(value / 0.0000152587890625);
    rgba.b = tmp / 255.0;
    value -= tmp * 0.0000152587890625;
    
    tmp = floor(value / 0.000000059604644775390625);
    rgba.g = tmp / 255.0;
    value -= tmp * 0.000000059604644775390625;
    
    tmp = floor(value / 0.00000000023283064365386962890625);
    rgba.r = tmp / 255.0;
    return rgba;
}

float RGBA2float(RGBA rgba){
    float tmp = 0.0;
    float sign = 1.0;
    if (rgba.a >= 128.0/255.0) {
        sign = -1.0;
        rgba.a -= 128.0/255.0;
    }
    tmp += rgba.a * 0.99609375;
    tmp += rgba.b * 0.0038909912109375;
    tmp += rgba.g * 0.000015199184417724609375;
    tmp += rgba.r * 0.00000005937181413173675537109375;
    //return tmp;
    return sign * (rgba.a * 256.0 * 256.0 * 256.0 * 255.0 + rgba.b * 256 * 256.0 * 255.0 + rgba.g * 256 * 255.0 + rgba.r * 255.0) / (256.0 * 256.0 * 256.0 * 128.0);
}

RGBA addRGBA(RGBA a, RGBA b){
    float tmp;
    float up = 0.00390625;
    float upword = 0.0;
    RGBA r;
    tmp = a.r + b.r;
    tmp *= 0.99609375;
    if(tmp > 0.99609375){
        r.r = tmp - 1.0;
        upword = up;
    }else{
        r.r = tmp;
    }
    tmp = a.g + b.g;
    tmp *= 0.99609375;
    tmp += upword;
    upword = 0.0;
    if(tmp > 0.99609375){
        r.g = tmp - 1.0;
        upword = up;
    }else{
        r.g = tmp;
    }
    tmp = a.b + b.b;
    tmp *= 0.99609375;
    tmp += upword;
    upword = 0.0;
    if(tmp > 0.99609375){
        r.b = tmp - 1.0;
        upword = up;
    }else{
        r.b = tmp;
    }
    tmp = a.a + b.a;
    tmp *= 0.99609375;
    tmp += upword;
    upword = 0.0;
    if(tmp > 0.99609375){
        r.a = tmp - 1.0;
        upword = up;
    }else{
        r.a = tmp;
    }
    
    r.a /= 0.99609375;
    r.b /= 0.99609375;
    r.g /= 0.99609375;
    r.r /= 0.99609375;
    return r;
}

void restrc(float* r1, int i1, int j1, float* r2, int i2, int j2){
    int i, j;
    for(int ii = 1; ii <= i2 - 1;ii++){
        i = 2 * ii;
        for (int jj = 1; jj <= j2 - 1; jj++) {
            j = 2 * jj;
            r2[jj * i2 + ii] = r1[j * i1 + i] + (r1[j * i1 + i - 1] +r1[j * i1 + i + 1] +r1[(j - 1) * i1 + i] +r1[(j + 1) * i1 + i]) / 2.0 + (r1[(j + 1) * i1 + i - 1] +r1[(j - 1) * i1 + i + 1] +r1[(j - 1) * i1 + i - 1] +r1[(j + 1) * i1 + i + 1]) / 4.0;
        }
    }
}

void prolon(float* d2, int i2, int j2, float* d1, int i1, int j1){
    for (int i = 1; i < i2 - 1; i++) {
        for (int j = 1; j < j2 - 1; j++) {
            d1[2 * j * i1 + 2 * i] = d2[j * i2 + i];
        }
    }
    int i, j;
    for (int ii = 1; ii < i2 - 1; ii++) {
        i = 2 * ii;
        d1[1 * i1 + i] = (3.0 * d1[0 * i1 + i] + 6.0 * d1[2 * i1 + i] - d1[4 * i1 + i]) / 8.0;
        for (j = 3; j < j1 - 3; j += 2) {
            d1[j * i1 + i] = (9.0 * d1[(j - 1) * i1 + i] + 9.0 * d1[(j + 1) * i1 + i] - d1[(j + 3) * i1 + i] - d1[(j - 3) * i1 + i]) / 16.0;
        }
        d1[(j1 - 2) * i1 + i] = (6.0 * d1[(j1 - 3) * i1 + i] + 3.0 * d1[(j1 - 1) * i1 + i] - d1[(j1 - 5) * i1 + i]);
    }
    for (int j = 1; j < j1 - 1; j++) {
        d1[j * i1 + 1] = (3.0 * d1[j * i1 + 0] + 6.0 * d1[j * i1 + 2] - d1[j * i1 + 4]) / 8.0;
        for (int i = 3; i < i1 - 3; i += 2) {
            d1[j * i1 + i] = (9.0 * d1[j * i1 + i - 1] + 9.0 * d1[j * i1 + i + 1] - d1[j * i1 + i + 3] - d1[j * i1 + i - 3]) / 16.0;
        }
        d1[j * i1 + i1 - 2] = (6.0 * d1[j * i1 + i1 - 3] + 3.0 * d1[j * i1 + i1 - 1] - d1[j * i1 + i1 - 5]);
    }
}

void prorlx(float* d2, int i2, int j2, float* d1, int i1, int j1, float* r1){
    for (int i = 1; i <= i2 - 1; i++) {
        for (int j = 1; j <= j2 - 1; j++) {
            d1[2 * j * i1 + 2 * i] = d2[j * i2 + i];
        }
    }
    int i;
    for (int ii = 1; ii <= i2 - 1; ii++) {
        i = 2 * ii;
        d1[1 * i1 + i] = (3.0 * d1[0 * i1 + i] + 6.0 * d1[2 * i1 + i] - d1[4 * i1 + i]) / 8.0;
        for (int j = 3; j <= j1 - 3; j += 2) {
            d1[j * i1 + i] = (9.0 * d1[(j - 1) * i1 + i] + 9.0 * d1[(j + 1) * i1 + i] - d1[(j + 3) * i1 + i] - d1[(j - 3) * i1 + i]) / 16.0;
        }
        d1[(j1 - 1) * i1 + i] =
        (6.0 * d1[(j1 - 2) * i1 + i]
         + 3.0 * d1[j1 * i1 + i]
         - d1[(j1 - 4) * i1 + i]);
    }
    for (int j = 1; j <= j1 - 1; j++) {
        d1[j * i1 + 1] = (3.0 * d1[j * i1 + 0] + 6.0 * d1[j * i1 + 2] - d1[j * i1 + 4]) / 8.0;
        for (int i = 3; i <= i1 - 3; i += 2) {
            d1[j * i1 + i] = (9.0 * d1[j * i1 + i - 1] + 9.0 * d1[j * i1 + i + 1] - d1[j * i1 + i + 3] - d1[j * i1 + i - 3]) / 16.0;
        }
        d1[j * i1 + i1 - 1] = (6.0 * d1[j * i1 + i1 - 2] + 3.0 * d1[j * i1 + i1] - d1[j * i1 + i1 - 4]);
    }
    for (int i = 1; i <= i1 - 1; i++) {
        for(int j = 1; j <= j1 - 1;j++){
            r1[j * i1 + i] = r1[j * i1 + i] + d1[j * i1 + i - 1] + d1[j * i1 + i + 1] + d1[(j + 1) * i1 + i] + d1[(j - 1) * i1 + i] - 4.0 * d1[j * i1 + i];
        }
    }
    for (int i = 1; i <= i1 - 1; i++) {
        for(int j = 1; j <= j1 - 1;j++){
            d1[j * i1 + i] = d1[j * i1 + i] + r1[j * i1 + i] / 4.0;
        }
    }
}

void gaussb(float* r1, float* d1, float* f1, int i1, int j1){
    int updated = 1;
    float precision = 0.0001;
    float prev_d1;
    float tmp;
    while (updated > 0) {
        updated = 0;
        for (int j = 1 ; j <= j1 - 1; j++)
        {
            for (int i = 1; i <= i1 - 1; i++)
            {
                prev_d1 = d1[j * i1 + i];
                tmp = f1[j * i1 + i];
                tmp = r1[j * i1 + i];
                d1[j * i1 + i] = 0.25 * (d1[j * i1 + i - 1] + d1[j * i1 + i + 1] + d1[(j - 1) * i1 + i] + d1[(j + 1) * i1 + i] + r1[j * i1 + i] - 4.0 * f1[j * i1 + i]);
                if (absd(d1[j * i1 + i] - prev_d1) > precision) {
                    updated++;
                }
            }
        }
        //NSLog(@"gaussb updated %d", updated);
    }
}

void zeros(float* p1, float* p2, int i1, int j1){
    for (int i = 0; i < i1; i++) {
        for (int j = 0; j < j1; j++) {
            p1[j * i1 + i] = 0.0;
            p2[j * i1 + i] = 0.0;
        }
    }
}

void solve(float* x1, float* f1, int i1, int j1){
    int updated = 1;
    float prev_i, current_err = 0.0, threshold = 0.001;
    
    while (updated > 0) {
        updated = 0;
        for (int j = 1 ; j < j1 - 1; j++)
        {
            for (int i = 1; i < i1 - 1; i++)
            {
                prev_i = x1[j * i1 + i];
                x1[j * i1 + i] = 0.25 * (x1[j * i1 + i + 1] + x1[j * i1 + i - 1] + x1[(j + 1) * i1 + i] + x1[(j - 1) * i1 + i] - f1[j * i1 + i]);
                //result[j * width + i] = prev_i + 1.25 * (result[j * width + i] - prev_i);

                current_err = absd(prev_i - x1[j * i1 + i]);
                if (current_err > threshold) {
                    updated++;
                }

            }
        }
        NSLog(@"Updated:%d", updated);
    }
    
}

void mgv(float* u1, float* f1, int i1, int j1, int current_grid, int max_grid){
    int i2 = i1 / 2;
    int j2 = j1 / 2;
    float* r1 = (float*)malloc(sizeof(float) * (i1 + 1) * (j1 + 1));
    float* d1 = (float*)malloc(sizeof(float) * (i1 + 1) * (j1 + 1));
    float* r2 = (float*)malloc(sizeof(float) * (i2 + 1) * (j2 + 1));
    float* d2 = (float*)malloc(sizeof(float) * (i2 + 1) * (j2 + 1));
    
    zeros(d1, r1, i1, j1);
    zeros(d2, r2, i2, j2);
    
    for (int i = 1; i < i1 - 1; i++) {
        for (int j = 1; j < j1 - 1; j++) {
            r1[j * i1 + i] = u1[j * i1 + i - 1] + u1[j * i1 + i + 1] + u1[(j - 1) * i1 + i] + u1[(j + 1) * i1 + i] - 4.0 * u1[j * i1 + i] - f1[j * i1 + i];
        }
    }

    restrc(r1, i1, j1, r2, i2, j2);
    
    if (current_grid == max_grid) {
        solve(d2, r2, i2, j2);
    }else{
        mgv(d2, r2, i2, j2, current_grid + 1, max_grid);
    }
    
    prolon(d2, i2, j2, d1, i1, j1);
    
    for (int i = 1; i < i1 - 1; i++) {
        for (int j = 1; j < j1 - 1; j++) {
            u1[j * i1 + i] += d1[j * i1 + i];
        }
    }
    
    free(d1);
    free(d2);
    free(r1);
    free(r2);
}

void fmg(float* u1, float* f1, int width, int height, int grids){
    int i, j;
    int i1 = width;
    int j1 = height;
    int i2 = i1 / 2;
    int j2 = j1 / 2;
    int i3 = i2 / 2;
    int j3 = j2 / 2;
    int i4 = i3 / 2;
    int j4 = j3 / 2;
    
    float* f2 = (float*)malloc(sizeof(float) * (i2 + 1) * (j2 + 1));
    float* u2 = (float*)malloc(sizeof(float) * (i2 + 1) * (j2 + 1));
    float* f3 = (float*)malloc(sizeof(float) * (i3 + 1) * (j3 + 1));
    float* u3 = (float*)malloc(sizeof(float) * (i3 + 1) * (j3 + 1));
    float* f4 = (float*)malloc(sizeof(float) * (i4 + 1) * (j4 + 1));
    float* u4 = (float*)malloc(sizeof(float) * (i4 + 1) * (j4 + 1));
    
    zeros(f2, u2, i2, j2);
    zeros(f3, u3, i3, j3);
    zeros(f4, u4, i4, j4);
    
    if (grids == 3) {
        restrc(f1, i1, j1, f2, i2, j2);
        restrc(f2, i2, j2, f3, i3, j3);
        restrc(f3, i3, j3, f4, i4, j4);
        solve(u4, f4, i4, j4);
        prolon(u4, i4, j4, u3, i3, j3);
        mgv(u3, f3, i3, j3, 1, 3);
        prolon(u3, i3, j3, u2, i2, j2);
        mgv(u2, f2, i2, j2, 2, 3);
        prolon(u2, i2, j2, u1, i1, j1);
        mgv(u1, f1, i1, j1, 3, 3);
    } else if (grids == 2){
        restrc(f1, i1, j1, f2, i2, j2);
        restrc(f2, i2, j2, f3, i3, j3);
        solve(u3, f3, i3, j3);
        prolon(u3, i3, j3, u2, i2, j2);
        mgv(u2, f2, i2, j2, 2, 3);
        prolon(u2, i2, j2, u1, i1, j1);
        mgv(u1, f1, i1, j1, 1, 3);
        solve(u1, f1, i1, j1);
    } else if (grids == 1){
        restrc(f1, i1, j1, f2, i2, j2);
        solve(u2, f2, i2, j2);
        prolon(u2, i2, j2, u1, i1, j1);
        mgv(u1, f1, i1, j1, 1, 2);
        solve(u1, f1, i1, j1);
    } else if(grids == 0){
        solve(u1, f1, i1, j1);
    }
    
    
    free(f2);
    free(u2);
    free(f3);
    free(u3);
    free(f4);
    free(u4);
}

void mgm2(float* u1, float* radiances, float* f1, int width, int height){
    int i, j;
    int i1 = width - 1;
    int j1 = height - 1;
    int i2 = i1 / 2;
    int j2 = j1 / 2;
    int i3 = i2 / 2;
    int j3 = j2 / 2;
    int i4 = i3 / 2;
    int j4 = j3 / 2;
    
    float rmean = 1.0;
    float tmean = 0.0001;
    
    float* r1 = (float*)malloc(sizeof(float) * (i1 + 1) * (j1 + 1));
    float* r2 = (float*)malloc(sizeof(float) * (i2 + 1) * (j2 + 1));
    float* r3 = (float*)malloc(sizeof(float) * (i3 + 1) * (j3 + 1));
    float* r4 = (float*)malloc(sizeof(float) * (i4 + 1) * (j4 + 1));
    float* d1 = (float*)malloc(sizeof(float) * (i1 + 1) * (j1 + 1));
    float* d2 = (float*)malloc(sizeof(float) * (i2 + 1) * (j2 + 1));
    float* d3 = (float*)malloc(sizeof(float) * (i3 + 1) * (j3 + 1));
    float* d4 = (float*)malloc(sizeof(float) * (i4 + 1) * (j4 + 1));
    float* f4 = (float*)malloc(sizeof(float) * (i4 + 1) * (j4 + 1));
    
    for (int ii = 0; ii <= i4 - 1; ii++) {
        for (int jj = 0; jj <= j4 - 1; jj++) {
            i = 4 * ii;
            j = 4 * jj;
            f4[jj * i4 + ii] = f1[j * i1 + i];
            float tmp = f1[j * i1 + i];
        }
    }
    
    float tmp = f1[24 * i1 + 0];
    tmp = f4[6 * i4 + 0];
    
    while (rmean > tmean) {
        rmean = 0.0;
        
        for (int i = 1; i <= i1 - 1; i++) {
            for (int j = 1; j <= j1 - 1; j++) {
                r1[j * i1 + i] = u1[j * i1 + i - 1] + u1[j * i1 + i + 1] + u1[(j - 1) * i1 + i] + u1[(j + 1) * i1 + i] - 4.0 * u1[j * i1 + i] - f1[j * i1 + i];
            }
        }
        restrc(r1, i1, j1, r2, i2, j2);
        restrc(r2, i2, j2, r3, i3, j3);
        restrc(r3, i3, j3, r4, i4, j4);
        gaussb(r4, d4, f4, i4, j4);
        prorlx(d4, i4, j4, d3, i3, j3, r3);
        prorlx(d3, i3, j3, d2, i2, j2, r2);
        prorlx(d2, i2, j2, d1, i1, j1, r1);
        
        for (int i = 1; i <= i1 - 1; i++) {
            for (int j = 1; j <= j1 - 1; j++) {
                u1[j * i1 + i] += d1[j * i1 + i];
                rmean += absd(r1[j * i1 + i]);
            }
        }
        rmean /= i1 * j1;
        NSLog(@"rmean: %lf", rmean);
        if(rmean > 200){
            break;
        }
        
    }
    
    free(r1);
    free(r2);
    free(r3);
    free(r4);
    free(d1);
    free(d2);
    free(d3);
    free(d4);
    free(f4);
}

- (UIImage *)hdrWithInputImage:(UIImage *)inputImage
{
    UIImage* resultImage;
    
    /*
    {
        float value = -2.3456789876543212345678987654321;
        value /= 10.0;
        RGBA rgba = float2RGBA(value);
        float result = RGBA2float(rgba);
        rgba = float2RGBA(result);
        result = RGBA2float(rgba);
        NSLog(@"input %lf, rgba:%lf %lf %lf %lf", value, rgba.a, rgba.b, rgba.g, rgba.r);
        NSLog(@"%lf", result);
        float err = result - value;
        NSLog(@"Err10:%lf", err * 100000000);
    }
    return inputImage;
    for (int i = 0; i < 10000; i++) {
        
        float value = (float)i / 10000;
        float unko = value / 2.0;
        RGBA v = float2RGBA(value);
        RGBA u = float2RGBA(unko);
        RGBA r = addRGBA(v, u);
        float result = RGBA2float(r);
        float tmp = value + unko;
        float err = result - tmp;

        NSLog(@"%lf -> %lf", tmp, result);
     
    }
    {
        float value = 9.9999999976;
        value /= 10.0;
        RGBA rgba = float2RGBA(value);
        float result = RGBA2float(rgba);
        rgba = float2RGBA(result);
        result = RGBA2float(rgba);
        NSLog(@"input %lf, rgba:%lf %lf %lf %lf", value, rgba.a, rgba.b, rgba.g, rgba.r);
        NSLog(@"%lf", result);
        float err = result - value;
        NSLog(@"Err10:%lf", err * 1000000000);
    }
    return inputImage;

    
     */
    
    
    int width = (int)CGImageGetWidth(inputImage.CGImage);
    int height = (int)CGImageGetHeight(inputImage.CGImage);
    size_t bitsPerComponent = CGImageGetBitsPerComponent(inputImage.CGImage);
    size_t bitsPerPixel = CGImageGetBitsPerPixel(inputImage.CGImage);
    size_t bytesPerRow = CGImageGetBytesPerRow(inputImage.CGImage);
    CGColorSpaceRef colorSpace = CGImageGetColorSpace(inputImage.CGImage);
    CGBitmapInfo bitmapInfo = CGImageGetBitmapInfo(inputImage.CGImage);
    BOOL shouldInterpolate = CGImageGetShouldInterpolate(inputImage.CGImage);
    CGColorRenderingIntent intent = CGImageGetRenderingIntent(inputImage.CGImage);
    
    // データプロバイダを取得する
    CGDataProviderRef dataProvider = CGImageGetDataProvider(inputImage.CGImage);
    
    // ビットマップデータを取得する
    CFDataRef data = CGDataProviderCopyData(dataProvider);
    CFMutableDataRef mutableData = CFDataCreateMutableCopy(0, 0, data);
    
    UInt8* buffer = (UInt8*)CFDataGetMutableBytePtr(mutableData);
    int i, j;
    UInt8* pixel;
    NSDate* start = [NSDate date];
    
    float r, g, b, _r, _g, _b, ev0, ev2, ev4, ev_2, ev_4, ev_8, exp0, exp2, exp4, exp8, exp_2, exp_4;
    
    exp0 = 1.0;
    exp2 = 0.5;
    exp4 = 0.25;
    exp_2 = 2.0;
    exp_4 = 4.0;
    
    float lum, _lum, alpha;
    float tmp;
    float *result = (float*)malloc(sizeof(float) * width * height);
    float *radiances = (float*)malloc(sizeof(float) * width * height * 4);
    float *g_div = (float*)malloc(sizeof(float) * width * height);
    float* phis = (float*)malloc(sizeof(float) * width * height);
    float *radiances_pixel;
    float *lum_delta;
    int radiance_index;
    float max_lum = 0.0;
    
    float max_pixel_value = (1.0 + 1.0 / exp2 + 1.0 / exp4 + 1.0 / exp_2 + 1.0 / exp_4) / 3.0;
    float min_pixel_value = 0.0;
    float lw;
    float lw_global;
    float l_white = 0.0;
    float lw_sum = 0.0;
    float w0, w1, w2, w3, w4;
    float set[3] = {0.0, 0.0, 0.0};
    
    YUV yuv;
    RGB rgb;
    
    zeros(result, g_div, width, height);

    
    for (j = 0 ; j < height; j++)
    {
        for (i = 0; i < width; i++)
        {
            // ピクセルのポインタを取得する
            pixel = buffer + j * bytesPerRow + i * 4;
            radiance_index = (j * width) * 4 + i * 4;
            //lum = pixel2luminance(pixel) / 255;
            
            
            r = (float)*(pixel + 0) / 255.0;
            g = (float)*(pixel + 1) / 255.0;
            b = (float)*(pixel + 2) / 255.0;
            
            set[0] = r;
            set[1] = g;
            set[2] = b;
            
            for(int i = 0;i < 3;i++){
                ev0 = set[i];
                ev2 = ev0 * ev0;
                ev4 = ev2 * ev2;
                ev_2 = ev0;
                screen(&ev_2);
                ev_4 = ev_2;
                screen(&ev_4);
                
                ev2 /= exp2;
                ev4 /= exp4;
                ev_2 /= exp_2;
                ev_4 /= exp_4;
                //w0 = ev_4 - ev4;
                //tmp = r + w0 * (ev2 + ev4 + ev_2 + ev_4) / 4.0;
                tmp = (ev0 + ev2 + ev4 + ev_2 + ev_4) / 5.0;
                //radiances[radiance_index + i] = tmp;
                
                // NUWAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
                radiances[radiance_index + i] = ev0;
            }
            
            lum = rgb2luminance(radiances[radiance_index], radiances[radiance_index + 1], radiances[radiance_index + 2]);
            lw_sum += log(0.0001 + lum);
            if(lum > l_white){
                l_white = lum;
            }
            radiances[radiance_index + 3] = log(lum + 0.0001);
        }
    }
    
    NSLog(@"called phis.");
    calcphis(phis, radiances, width, height, 0.1, 0.8);
    NSLog(@"called gradiant.");
    gradient(g_div, radiances, phis, width, height);
    free(phis);
    NSLog(@"called gassb.");
    //gaussb_(result, g_div, width, height);
    fmg(result, g_div, width, height, 2);
    free(g_div);
    NSLog(@"called expall.");
    expall(result, radiances, width, height);

    
    float s = 0.6;
    alpha = 0.3;
    
    
    for (int j = 1 ; j < height; j++)
    {
        for (int i = 1; i < width; i++)
        {
            radiance_index = (j * width) * 4 + i * 4;
            pixel = buffer + j * bytesPerRow + i * 4;
            
            tmp = result[j * width + i] / 1.0;
            
            
            rgb.r = pow(radiances[radiance_index + 0] / radiances[radiance_index + 3], s) * tmp;
            rgb.g = pow(radiances[radiance_index + 1] / radiances[radiance_index + 3], s) * tmp;
            rgb.b = pow(radiances[radiance_index + 2] / radiances[radiance_index + 3], s) * tmp;
            
            
            /*
            
            // Soft Light
            tmp = (float)*(pixel) / 255.0;
            tmp = rgb.r * alpha + (1.0 - alpha) * tmp;
            if(rgb.r < 0.5){
                rgb.r = pow(tmp, (1.0 - rgb.r) / 0.5);
            }else{
                rgb.r = pow(tmp, 0.5 / rgb.r);
            }
            tmp = (float)*(pixel + 1) / 255.0;
            tmp = rgb.g * alpha + (1.0 - alpha) * tmp;
            if(rgb.g < 0.5){
                rgb.g = pow(tmp, (1.0 - rgb.g) / 0.5);
            }else{
                rgb.g = pow(tmp, 0.5 / rgb.g);
            }
            tmp = (float)*(pixel + 2) / 255.0;
            tmp = rgb.b * alpha + (1.0 - alpha) * tmp;
            if(rgb.b < 0.5){
                rgb.b = pow(tmp, (1.0 - rgb.b) / 0.5);
            }else{
                rgb.b = pow(tmp, 0.5 / rgb.b);
            }
             
            
            
            
             // Hard light
             tmp = (float)*(pixel) / 255.0;
             if(rgb.r < 0.5){
             rgb.r = tmp * rgb.r * 2.0;
             }else{
             rgb.r = 2.0 * (rgb.r + tmp - tmp * rgb.r) - 1.0;
             }
             tmp = (float)*(pixel + 1) / 255.0;
             if(rgb.g < 0.5){
             rgb.g = tmp * rgb.g * 2.0;
             }else{
             rgb.g = 2.0 * (rgb.g + tmp - tmp * rgb.g) - 1.0;
             }
             tmp = (float)*(pixel + 2) / 255.0;
             if(rgb.b < 0.5){
             rgb.b = tmp * rgb.b * 2.0;
             }else{
             rgb.b = 2.0 * (rgb.b + tmp - tmp * rgb.b) - 1.0;
             }
             */
            
            
            
            
            //rgb.r = rgb.g = rgb.b = exp(gausian_pyramids[width * height * 4 + j * width + i]);
            //rgb.r = rgb.g = rgb.b = result[j * width + i];
            
            *(pixel) = (int)MAX(MIN(round(rgb.r * 255.0), 255.0), 0.0);
            *(pixel + 1) = (int)MAX(MIN(round(rgb.g * 255.0), 255.0), 0.0);
            *(pixel + 2) = (int)MAX(MIN(round(rgb.b * 255.0), 255.0), 0.0);
        }
    }

    free(result);
    free(radiances);
    
    // do stuff...
    NSTimeInterval timeInterval = [start timeIntervalSinceNow];
    NSLog(@"time:%lf", -timeInterval);
    
    
    // 効果を与えたデータを作成する
    CFDataRef effectedData;
    effectedData = CFDataCreate(NULL, buffer, CFDataGetLength(mutableData));
    CFRelease(mutableData);
    CFRelease(data);
    
    // 効果を与えたデータプロバイダを作成する
    CGDataProviderRef effectedDataProvider;
    effectedDataProvider = CGDataProviderCreateWithCFData(effectedData);
    
    // 画像を作成する
    CGImageRef effectedCgImage = CGImageCreate(
                                               width, height,
                                               bitsPerComponent, bitsPerPixel, bytesPerRow,
                                               colorSpace, bitmapInfo, effectedDataProvider,
                                               NULL, shouldInterpolate, intent);
    
    
    resultImage = [[UIImage alloc] initWithCGImage:effectedCgImage];
    
    // 作成したデータを解放する
    CGImageRelease(effectedCgImage);
    CFRelease(effectedDataProvider);
    CFRelease(effectedData);
    return resultImage;
}


- (UIImage *)process:(UIImage *)inputImage
{
    UIImage* resultImage;
    
    int width = (int)CGImageGetWidth(inputImage.CGImage);
    int height = (int)CGImageGetHeight(inputImage.CGImage);
    size_t bitsPerComponent = CGImageGetBitsPerComponent(inputImage.CGImage);
    size_t bitsPerPixel = CGImageGetBitsPerPixel(inputImage.CGImage);
    size_t bytesPerRow = CGImageGetBytesPerRow(inputImage.CGImage);
    CGColorSpaceRef colorSpace = CGImageGetColorSpace(inputImage.CGImage);
    CGBitmapInfo bitmapInfo = CGImageGetBitmapInfo(inputImage.CGImage);
    BOOL shouldInterpolate = CGImageGetShouldInterpolate(inputImage.CGImage);
    CGColorRenderingIntent intent = CGImageGetRenderingIntent(inputImage.CGImage);
    
    // データプロバイダを取得する
    CGDataProviderRef dataProvider = CGImageGetDataProvider(inputImage.CGImage);
    
    // ビットマップデータを取得する
    CFDataRef data = CGDataProviderCopyData(dataProvider);
    CFMutableDataRef mutableData = CFDataCreateMutableCopy(0, 0, data);
    
    UInt8* buffer = (UInt8*)CFDataGetMutableBytePtr(mutableData);
    int i, j;
    UInt8* pixel;
    NSDate* start = [NSDate date];
    
    float r, g, b, _r, _g, _b, ev0, ev2, ev4, ev8, ev_2, ev_4, ev_8, exp0, exp2, exp4, exp8, exp_2, exp_4, exp_8;
    
    exp0 = 1.0;
    exp2 = 0.5;
    exp4 = 0.25;
    exp8 = 0.125;
    exp_2 = 2.0;
    exp_4 = 4.0;
    exp_8 = 8.0;
    
    float lum, _lum;
    float tmp;
    float *radiances;
    float *radiances_pixel;
    float *lum_delta;
    int radiance_index;
    
    float max_pixel_value = (1.0 + 1.0 / exp2 + 1.0 / exp4 + 1.0 / exp_2 + 1.0 / exp_4) / 3.0;
    float min_pixel_value = 0.0;
    float lw;
    float lw_global;
    float l_white = 0.0;
    float lw_sum = 0.0;
    float w0, w1, w2, w3, w4;
    float set[3] = {0,0, 0.0, 0.0};
    
    
    radiances = (float*)malloc(sizeof(float) * width * height * 4);
    lum_delta = (float*)malloc(sizeof(float) * width * height);

    for (j = 0 ; j < height; j++)
    {
        for (i = 0; i < width; i++)
        {
            // ピクセルのポインタを取得する
            pixel = buffer + j * bytesPerRow + i * 4;
            radiance_index = (j * width) * 4 + i * 4;
            //lum = pixel2luminance(pixel) / 255;
            
            
            r = (float)*(pixel + 0) / 255.0;
            g = (float)*(pixel + 1) / 255.0;
            b = (float)*(pixel + 2) / 255.0;
            
            set[0] = r;
            set[1] = g;
            set[2] = b;
            
            for(int i = 0;i < 3;i++){
                ev0 = set[i];
                ev2 = ev0 * ev0;
                ev4 = ev2 * ev2;
                ev_2 = ev0;
                screen(&ev_2);
                ev_4 = ev_2;
                screen(&ev_4);
                
                ev2 /= exp2;
                ev4 /= exp4;
                ev_2 /= exp_2;
                ev_4 /= exp_4;
                //w0 = ev_4 - ev4;
                //tmp = r + w0 * (ev2 + ev4 + ev_2 + ev_4) / 4.0;
                tmp = (ev0 + ev2 + ev4 + ev_2 + ev_4) / 5.0;
                radiances[radiance_index + i] = tmp;
            }
            
            lum = rgb2luminance(radiances[radiance_index], radiances[radiance_index + 1], radiances[radiance_index + 2]);
            lw_sum += log(0.0001 + lum);
            if(lum > l_white){
                l_white = lum;
            }
            radiances[radiance_index + 3] = log(lum + 0.0001);
        }
    }
    
    /*
    for (j = 0 ; j < height; j++)
    {
        for (i = 0; i < width; i++)
        {
            // ピクセルのポインタを取得する
            pixel = buffer + j * bytesPerRow + i * 4;
            radiance_index = (j * width) * 4 + i * 4;
            lum = pixel2luminance(pixel) / 255;
            
            
            r = (float)*(pixel + 0) / 255.0;
            g = (float)*(pixel + 1) / 255.0;
            b = (float)*(pixel + 2) / 255.0;
            if(lum > 0.5){
                ev2 = r * r;
                ev4 = ev2 * ev2;
                tmp = r + ev2 / exp2 + ev4 / exp4;
                radiances[radiance_index] = tmp / 3.0;
                ev2 = g * g;
                ev4 = ev2 * ev2;
                tmp = g + ev2 / exp2 + ev4 / exp4;
                radiances[radiance_index + 1] = tmp / 3.0;
                ev2 = b * b;
                ev4 = ev2 * ev2;
                tmp = b + ev2 / exp2 + ev4 / exp4;
                radiances[radiance_index + 2] = tmp / 3.0;
            }else{
                ev_2 = r;
                screen(&ev_2);
                ev_4 = ev_2;
                screen(&ev_4);
                tmp = r + ev_2 / exp_2 + ev_4 / exp_4;
                radiances[radiance_index] = tmp / 3.0;
                ev_2 = g;
                screen(&ev_2);
                ev_4 = ev_2;
                screen(&ev_4);
                tmp = g + ev_2 / exp_2 + ev_4 / exp_4;
                radiances[radiance_index + 1] = tmp / 3.0;
                ev_2 = b;
                screen(&ev_2);
                ev_4 = ev_2;
                screen(&ev_4);
                tmp = b + ev_2 / exp_2 + ev_4 / exp_4;
                radiances[radiance_index + 2] = tmp / 3.0;
            }
            lum = rgb2luminance(radiances[radiance_index], radiances[radiance_index + 1], radiances[radiance_index + 2]);
            lw_sum += lum;
            if(lum > l_white){
                l_white = lum;
            }
            radiances[radiance_index + 3] = lum;
        }
    }
*/
    
    lw_global = lw_sum / (width * height);
    lw_global = exp(lw_global);
    NSLog(@"Lw: %lf", lw);
    NSLog(@"L_white: %lf", l_white);
    NSLog(@"Lw_global: %lf", lw_global);
    float alpha = 0.18;
    float ratio = 1.0 / l_white;
    
    float l = 1.0;
    int pyramid_count = 0;
    int base_size = (width > height) ? width : height;
    while (base_size > 64) {
        base_size /= 2;
        pyramid_count++;
    }
    float* gausian_pyramids = (float *)malloc(sizeof(float) * width * height * pyramid_count);
    
    /*
    int cf = (int)pow(2.0, pyramid_count - 1) * 2 + 1;
    float* coefficient = (float*)malloc(sizeof(float) * (cf * (cf + 1)) / 2);
    coefficient[0] = 1.0;
    coefficient[1] = 1.0;
    coefficient[2] = 1.0;
    int index = 0;
    float sum;
    for(int i = 2;i < cf;i++){
        index = i * (i + 1) / 2;
        for(int j = 0;j <= i;j++){
            if(j == 0){
                coefficient[index + j] = 1.0;
            }else if(j == i){
                coefficient[index + j] = 1.0;
            }else{
                coefficient[index + j] = coefficient[index + j - i] + coefficient[index + j - i - 1];
            }
        }
    }
    for(int i = 2;i < cf;i++){
        index = i * (i + 1) / 2;
        sum = 0.0;
        for(int j = 0;j <= i;j++){
            sum += coefficient[index + j];
        }
        for(int j = 0;j <= i;j++){
            coefficient[index + j] /= sum;
        }
    }

    int range, x, y;
    
    for(int n = 0; n < pyramid_count;n++){
        range = (int)pow(2.0, n);
        index = (range * 2 + 1) * (range * 2 + 2) / 2 - range - 1;
        NSLog(@"middle:%d", index);
        for (int j = 1 ; j < height - 1; j++)
        {
            for (int i = 1; i < width - 1; i++)
            {
                sum = 0.0;
                for(int xx = -range;xx <= range;xx++){
                    x = i + xx;
                    if(x < 0){
                        x = -x;
                    }else if(x > width){
                        x = i - xx;
                    }
                    tmp = coefficient[index + xx];
                    sum += radiances[j * width * 4 + x * 4 + 3] * tmp;
                }
                for(int yy = -range;yy <= range;yy++){
                    y = j + yy;
                    if(y < 0){
                        y = -y;
                    }else if(y > height){
                        y = j - yy;
                    }
                    tmp = coefficient[index + yy];
                    sum += radiances[y * width * 4 + i * 4 + 3] * tmp;
                }
                gausian_pyramids[j * width * pyramid_count + i * pyramid_count + n] = sum;
            }
        }

    }

     */
    
    NSLog(@"pyramid_count:%d", pyramid_count);
    int range = 2;
    int _range = range;
    float N;
    float weight, sum;
    int x, y;
    float* g_table = (float*)malloc(sizeof(float) * width * height);
    float sigma;
    int w_table_length = (int)(range * 2.0 * (pow(2.0, (float)pyramid_count) - 1.0)) + pyramid_count;
    int w_table_index = 0;
    float* w_table = (float*)malloc(sizeof(float) * w_table_length);
    for(int i = 0;i < w_table_length;i++){
        w_table[i] = 77.7;
    }
    
    for(int n = 0;n < pyramid_count;n++){
        sigma = 0.3 * ((float)(_range) / 2.0 - 1.0) + 0.8;
        //sigma = range;
        for (y = 0; y < height; y++) {
            for(x = 0;x < width;x++){
                N = 0.0;
                sum = 0.0;
                for(int xx = -_range;xx <= _range;xx++){
                    if(n == 3){
                        ;
                    }
                    if(x + xx < 0 || x + xx >= width){
                        continue;
                    }
                    _lum = radiances[y * width * 4 + (x + xx) * 4 + 3];
                    w_table_index = (int)(range * 2.0 * (pow(2.0, (float)n + 1) - 1.0)) + n - _range + xx;
                    weight = w_table[w_table_index];
                    if(weight == 77.7){
                        weight = exp(-1 * ((xx * xx) / (2.0 * sigma * sigma)));
                        w_table[w_table_index] = weight;
                    }
                    N += weight;
                    sum += _lum * weight;
                }
                g_table[y * width + x] = sum / N;
            }
        }
        
        for(x = 0;x < width;x++){
            for (y = 0; y < height; y++) {
                N = 0.0;
                sum = 0.0;
                for(int yy = -_range;yy <= _range;yy++){
                    if(y + yy < 0 || y + yy >= height){
                        continue;
                    }
                    _lum = g_table[(y + yy) * width + x];
                    w_table_index = (int)(range * 2.0 * (pow(2.0, (float)n + 1) - 1.0)) + n - _range + yy;
                    weight = w_table[w_table_index];
                    if(weight == 77.7){
                        weight = exp(-1 * ((yy * yy) / (2.0 * sigma * sigma)));
                        w_table[w_table_index] = weight;
                    }
                    N += weight;
                    sum += _lum * weight;
                }
                gausian_pyramids[width * height * n + y * width + x] = sum / N;
            }
        }
        _range *= 2;
        NSLog(@"pyramid");
    }
    
    free(g_table);
    free(w_table);
    
    
    
    YUV yuv;
    RGB rgb;
    

    float* h_nabla = (float*)malloc(sizeof(float) * width * height * 2);
    float* g_div = (float*)malloc(sizeof(float) * width * height);
    float* phis = (float*)malloc(sizeof(float) * width * height);
    float* result = (float*)malloc(sizeof(float) * width * height);
    float h_x;
    float h_y;
    float h;
    float phi;
    float phi_max = 0.0;
    float L = 1.0;
    
    float g_div_avg = 0.0;
    
    
    for (int j = 0 ; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            h_nabla[j * width * 2 + i * 2 + 0] = 0.0;
            h_nabla[j * width * 2 + i * 2 + 1] = 0.0;
            g_div[j * width + i] = 0.0;
            if(j == 0 || j == height - 1){
                continue;
            }
            if(i == 0 || i == width - 1){
                continue;
            }
            
            // ピクセルのポインタを取得する
            radiance_index = (j * width) * 4 + i * 4;
            lum = radiances[radiance_index + 3];
            
            h_nabla[j * width * 2 + i * 2 + 0] = radiances[(j * width) * 4 + (i + 1) * 4 + 3] - lum;    //nabla Hx
            h_nabla[j * width * 2 + i * 2 + 1] = radiances[((j + 1) * width) * 4 + i * 4 + 3] - lum;    //nabla Hy
            tmp = sqrt(h_nabla[j * width * 2 + i * 2 + 0] * h_nabla[j * width * 2 + i * 2 + 0] + h_nabla[j * width * 2 + i * 2 + 1] * h_nabla[j * width * 2 + i * 2 + 1]);
            g_div_avg += tmp;
        }
    }
    
    g_div_avg /= width * height;
    NSLog(@"g_div_avg:%lf", g_div_avg);

    alpha = 0.1 * g_div_avg;
    float beta = 0.8;

    
    for (int j = 1 ; j < height - 1; j++)
    {
        for (int i = 1; i < width - 1; i++)
        {
            phi = 1.0;
            for(int n = pyramid_count - 1;n >= 0;n--){
                h_x = (gausian_pyramids[width * height * n + j * width + (i + 1)] - gausian_pyramids[width * height * n + j * width + (i - 1)]) / pow(2.0, (float)(n + 2));
                h_y = (gausian_pyramids[width * height * n + (j + 1) * width + i] - gausian_pyramids[width * height * n + (j - 1) * width + i]) / pow(2.0, (float)(n + 2));
                h = sqrt(h_x * h_x + h_y * h_y);
                tmp = L * (alpha / h) * pow(h / alpha, beta);
                if(tmp > 1.0){
                    tmp = 1.0;
                }
                phi *= tmp;
            }
            h_x = (radiances[j * width * 4 + (i + 1) * 4 + 3] - radiances[j * width * 4 + (i - 1) * 4 + 3]) / 2.0;
            h_y = (radiances[(j + 1) * width * 4 + i * 4 + 3] - radiances[(j - 1) * width * 4 + i * 4 + 3]) / 2.0;
            h = sqrt(h_x * h_x + h_y * h_y);
            
            tmp = (alpha / h) * pow(h / alpha, beta);
            if(tmp > 1.0){
                tmp = 1.0;
            }
            phi *= tmp;

            phis[j * width + i] = 0.0;
            if(isnan(phi)){
                phi = 1.0;
            }
            if(phi < 0.85){
                phi = 0.85;
            }
            phis[j * width + i] = phi;
        }
    }
        
    bool* flags = (bool*)malloc(sizeof(bool) * width * height);
    float max_phi = 1.0;
    float min_phi = 1.0;
    for (int j = 0 ; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {
            flags[j * width + i] = true;
            if(j == 0 || j == height - 1){
                continue;
            }
            if(i == 0 || i == width - 1){
                continue;
            }
            
            radiance_index = (j * width) * 4 + i * 4;
            lum = radiances[radiance_index + 3];
            
            
            phi = phis[j * width + i];
            if(max_phi < phi){
                max_phi = phi;
            }
            if(min_phi > phi){
                min_phi = phi;
            }
            //phi = 1.0;
            //h_nabla[j * width + i + 0] *= phi;
            //h_nabla[j * width + i + 1] *= phi;
            
            radiance_index = (j * width) * 4 + i * 4;
            lum = radiances[radiance_index + 3];
            
            g_div[j * width + i] += phi * (radiances[j * width * 4 + (i - 1) * 4 + 3] + radiances[j * width * 4 + (i + 1) * 4 + 3] - 2.0 * lum);
            g_div[j * width + i] += phi * (radiances[(j - 1) * width * 4 + i * 4 + 3] + radiances[(j + 1) * width * 4 + i * 4 + 3] - 2.0 * lum);
            
            tmp = g_div[j * width + i];
            /*
             g_div[j * width + i] += h_nabla[j * width * 2 + i * 2 + 0] - h_nabla[j * width * 2 + (i - 1) * 2 + 0];
             g_div[j * width + i] += h_nabla[j * width * 2 + i * 2 + 1] - h_nabla[(j - 1) * width * 2 + i * 2 + 1];
             */
            
            result[j * width + i] = 0.0;
        }
    }
    
    float prev_i;
    float threshold = 0.1;
    float max_i = 0.001;
    float current_err = 0.0;
    float max_err = 0.0;
    int updated = 1;
    bool cflag = true;
    int I, J;
    
    int widthDiv2 = ceil((float)width / 2.0);
    int heightDiv2 = ceil((float)height / 2.0);
    
    float* D = malloc(sizeof(float) * width * height);
    float* D1 = malloc(sizeof(float) * widthDiv2 * heightDiv2);
    float* R1 = malloc(sizeof(float) * widthDiv2 * heightDiv2);
    
    
    while (threshold > 0.001) {
        threshold /= 10.0;
        NSLog(@"threshold:%lf", threshold);
        updated = 1;
        I = 0;
        J = 0;
        for (int j = 2 ; j < height - 2; j += 2)
        {
            for (int i = 2; i < width - 2; i += 2)
            {
                R1[J * widthDiv2 + I] = result[j * width + i] + (result[j * width + i - 1] + result[j * width + i + 1] + result[(j - 1) * width + i] + result[(j + 1) * width + i]) / 2.0 + (result[(j + 1) * width + i - 1] + result[(j - 1) * width + i] + result[(j - 1) * width - i] + result[(j + 1) * width + i]) / 4.0;
                D1[J * widthDiv2 + I] = 0.0;
                D[j * width + i] = 0.0;
                I++;
            }
            J++;
            I = 0;
        }
        while (updated > 0) {
            updated = 0;
            I = 1;
            J = 1;
            for (int j = 0 ; j < height; j += 2)
            {
                if(J >= heightDiv2 - 1){
                    continue;
                }
                for (int i = 0; i < width; i += 2)
                {
                    if(I >= widthDiv2 - 1){
                        continue;
                    }
                    prev_i = D1[J * widthDiv2 + I];
                    D1[J * widthDiv2 + I] = 0.25 * (R1[J * widthDiv2 + I] + D1[J * widthDiv2 + I - 1] + D1[J * widthDiv2 + I + 1] + D1[(J + 1) * widthDiv2 + I] + D1[(J - 1) * widthDiv2 + I] - g_div[j * width + i]);
                    D[j * width + i] = D1[J * widthDiv2 + I];
                    if(absd(D1[J * widthDiv2 + I]) > max_i){
                        max_i = absd(D1[J * widthDiv2 + I]);
                    }
                    current_err = absd(prev_i - D1[J * widthDiv2 + I]) / max_i;
                    if(current_err > threshold){
                        updated++;
                    }else{
                    }
                    I++;
                }
                J++;
                I = 1;
            }
            NSLog(@"D1:Updated:%d", updated);

        }
        int i, j;
        for(int ii = 1;ii < widthDiv2;ii++){
            i = 2 * ii;
            D[1 * width + i] = (3.0 * D[0 * width + i] + 6.0 * D[2 * width + i] - D[4 * width + i]) / 8.0;
            for (j = 3; j < height - 3; j += 2)
            {
                D[j * width + i] = (9.0 * D[(j - 1) * width + i] + 9.0 * D[(j + 1) * width + i] - D[(j - 3) * width + i] - D[(j + 3) * width + i]) / 16.0;
            }
            j--;
            if(j % 2 == 0){
                j += 1;
            }else{
                j += 2;
            }
            D[j * width + i] = (3.0 * D[(j + 1) * width + i] + 6.0 * D[(j - 1) * width + i] - D[(j - 3) * width + i]) / 8.0;
        }
        for (int j = 0 ; j < height; j++)
        {
            i = 3;
            D[j * width + 1] = (3.0 * D[j * width + 0] + 6.0 * D[j * width + 2] - D[j * width + 4]) / 8.0;
            for (i = 3; i < width - 3; i += 2)
            {
                D[j * width + i] = (9.0 * D[j * width + i - 1] + 9.0 * D[j * width + i + 1] - D[j * width + i - 3] - D[j * width + i + 3]) / 16.0;
            }
            i--;
            if(i % 2 == 0){
                i += 1;
            }else{
                i += 2;
            }
            D[j * width + i] = (3.0 * D[j * width + i + 1] + 6.0 * D[j * width + i - 1] - D[j * width + i - 3]) / 8.0;
        }

        updated = 1;
        max_i = 0.0;
        cflag = !cflag;
        while (updated > 0) {
            updated = 0;
            for (int j = 1 ; j < height - 1; j++)
            {
                for (int i = 1; i < width - 1; i++)
                {
                    if(flags[j * width + i] == cflag){
                        continue;
                    }
                    
                    prev_i = result[j * width + i];
                    result[j * width + i] = 0.25 * (result[j * width + i + 1] + result[j * width + i - 1] + result[(j + 1) * width + i] + result[(j - 1) * width + i] - g_div[j * width + i]);
                    //result[j * width + i] = prev_i + 1.25 * (result[j * width + i] - prev_i);
                    if(absd(result[j * width + i]) > max_i){
                        max_i = absd(result[j * width + i]);
                    }
                    //tmp = result[(j + 1) * width + i] + result[(j - 1) * width + i] + result[j * width + (i + 1)] + result[j * width + (i - 1)] - 4.0 * result[j * width + i];
                    current_err = absd(prev_i - result[j * width + i]) / max_i;
                    if(current_err > threshold){
                        updated++;
                        flags[j * width + i] = !cflag;
                        flags[j * width + i - 1] = !cflag;
                        flags[j * width + i + 1] = !cflag;
                        flags[(j + 1) * width + i] = !cflag;
                        flags[(j - 1) * width + i] = !cflag;
                        flags[(j + 1) * width + i + 1] = !cflag;
                        flags[(j + 1) * width + i - 1] = !cflag;
                        flags[(j - 1) * width + i + 1] = !cflag;
                        flags[(j - 1) * width + i - 1] = !cflag;
                    }else{
                        flags[j * width + i] = cflag;
                    }
                }
            }
            NSLog(@"Updated:%d", updated);
        }
        
        
        for (int j = 1 ; j < height - 1; j++)
        {
            for (int i = 1; i < width - 1; i++)
            {
                result[j * width + i] += D[j * width + i + 1] + D[j * width + i - 1] + D[(j + 1) * width + i] + D[(j - 1) * width + i] - 4.0 * D[j * width + i];
                D[j * width + i] += result[j * width + i] / 4.0;
            }
        }
        
    }
    
    free(R1);
    free(D1);
    free(D);
    
    NSLog(@"max_phi:%lf", max_phi);
    NSLog(@"min_phi:%lf", min_phi);

    
    float* new_h_nabla = (float*)malloc(sizeof(float) * width * height * 2);
    float* new_g_div = (float*)malloc(sizeof(float) * width * height);
    
    for (int j = 0 ; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {

            new_h_nabla[j * width * 2 + i * 2 + 0] = 0.0;
            new_h_nabla[j * width * 2 + i * 2 + 1] = 0.0;
            new_g_div[j * width + i] = 0.0;
            if(j == 0 || j == height - 1){
                continue;
            }
            if(i == 0 || i == width - 1){
                continue;
            }

            lum = result[j * width + i];
            
            new_g_div[j * width + i] += result[j * width + i + 1] + result[j * width + i - 1] - 2.0 * lum;
            new_g_div[j * width + i] += result[(j + 1) * width + i] + result[(j - 1) * width + i] - 2.0 * lum;
            
            
            if(i == 194 && j == 267){
                NSLog(@"original g_div:%lf", g_div[j * width + i]);
                NSLog(@"new g_div:%lf", new_g_div[j * width + i]);
                NSLog(@"original lum:%lf", radiances[j * width * 4 + i * 4 + 3]);
                NSLog(@"new lum:%lf", result[j * width + i]);
                NSLog(@"phi:%lf", phis[j * width + i]);
            }
            
        }
    }
    l_white = 0.0;
    float l_min = 1.0;
    
    for (int j = 1 ; j < height - 1; j++)
    {
        for (int i = 1; i < width - 1; i++)
        {
            radiance_index = (j * width) * 4 + i * 4;
            tmp = result[j * width + i];
            result[j * width + i] = exp(tmp);
            tmp = radiances[radiance_index + 3];
            radiances[radiance_index + 3] = exp(tmp);
            if(result[j * width + i] > l_white){
                l_white = result[j * width + i];
            }
            
            if(result[j * width + i] < l_min){
                l_min = result[j * width + i];
            }
        }
    }
    
    float diff;
    float max_diff = 0.0;
    float min_diff = 1.0;
    float max_lum = 0.0;
    
    for (int j = 1 ; j < height - 1; j++)
    {
        for (int i = 1; i < width - 1; i++)
        {
            
            radiance_index = (j * width) * 4 + i * 4;
            tmp = result[j * width + i] / (1.0 + result[j * width + i]);
            result[j * width + i] = tmp;
            lum = radiances[radiance_index + 3];
            
            diff = tmp - lum;
            if(diff > 0.0){
                diff *= 6.0 * lum / (1.0 + 6.0 * lum);
                //result[j * width + i] = lum + diff;
            }
                        
            if(diff > max_diff){
                max_diff = diff;
            }
            if(diff < min_diff){
                min_diff = diff;
            }
            if(tmp > max_lum){
                max_lum = tmp;
            }
            
            //result[j * width + i] = 4.0 * lum / (1.0 + 4.0 * lum);
        }
    }
    
    NSLog(@"max_diff:%lf", max_diff);
    NSLog(@"min_diff:%lf", min_diff);
    NSLog(@"max_lum:%lf", max_lum);
    
    NSLog(@"l_white:%lf", l_white);
    NSLog(@"l_min:%lf", l_min);
    
    ratio = 1.0 / l_white;
    float s = 0.6;
    alpha = 0.3;

    
    for (int j = 1 ; j < height - 1; j++)
    {
        for (int i = 1; i < width - 1; i++)
        {
            radiance_index = (j * width) * 4 + i * 4;
            pixel = buffer + j * bytesPerRow + i * 4;
            
            tmp = result[j * width + i] / max_lum;
    
            
            rgb.r = pow(radiances[radiance_index + 0] / radiances[radiance_index + 3], s) * tmp;
            rgb.g = pow(radiances[radiance_index + 1] / radiances[radiance_index + 3], s) * tmp;
            rgb.b = pow(radiances[radiance_index + 2] / radiances[radiance_index + 3], s) * tmp;
            
            
            
            // Soft Light
            tmp = (float)*(pixel) / 255.0;
            tmp = rgb.r * alpha + (1.0 - alpha) * tmp;
            if(rgb.r < 0.5){
                rgb.r = pow(tmp, (1.0 - rgb.r) / 0.5);
            }else{
                rgb.r = pow(tmp, 0.5 / rgb.r);
            }
            tmp = (float)*(pixel + 1) / 255.0;
            tmp = rgb.g * alpha + (1.0 - alpha) * tmp;
            if(rgb.g < 0.5){
                rgb.g = pow(tmp, (1.0 - rgb.g) / 0.5);
            }else{
                rgb.g = pow(tmp, 0.5 / rgb.g);
            }
            tmp = (float)*(pixel + 2) / 255.0;
            tmp = rgb.b * alpha + (1.0 - alpha) * tmp;
            if(rgb.b < 0.5){
                rgb.b = pow(tmp, (1.0 - rgb.b) / 0.5);
            }else{
                rgb.b = pow(tmp, 0.5 / rgb.b);
            }
            
            
            
            /*
            // Hard light
            tmp = (float)*(pixel) / 255.0;
            if(rgb.r < 0.5){
                rgb.r = tmp * rgb.r * 2.0;
            }else{
                rgb.r = 2.0 * (rgb.r + tmp - tmp * rgb.r) - 1.0;
            }
            tmp = (float)*(pixel + 1) / 255.0;
            if(rgb.g < 0.5){
                rgb.g = tmp * rgb.g * 2.0;
            }else{
                rgb.g = 2.0 * (rgb.g + tmp - tmp * rgb.g) - 1.0;
            }
            tmp = (float)*(pixel + 2) / 255.0;
            if(rgb.b < 0.5){
                rgb.b = tmp * rgb.b * 2.0;
            }else{
                rgb.b = 2.0 * (rgb.b + tmp - tmp * rgb.b) - 1.0;
            }
             */
            
             
             
            
            //rgb.r = rgb.g = rgb.b = exp(gausian_pyramids[width * height * 4 + j * width + i]);
            //rgb.r = rgb.g = rgb.b = result[j * width + i];
            
            *(pixel) = (int)MAX(MIN(round(rgb.r * 255.0), 255.0), 0.0);
            *(pixel + 1) = (int)MAX(MIN(round(rgb.g * 255.0), 255.0), 0.0);
            *(pixel + 2) = (int)MAX(MIN(round(rgb.b * 255.0), 255.0), 0.0);
        }
    }

    
    free(new_g_div);
    free(new_h_nabla);


    free(g_div);
    free(result);
    free(phis);
    free(h_nabla);
    free(gausian_pyramids);
    free(flags);
    //free(new_g_div);
    //free(new_h_nabla);
    
    
    /*
    
    for (j = 0 ; j < height; j++)
    {
        for (i = 0; i < width; i++)
        {
            // ピクセルのポインタを取得する
            pixel = buffer + j * bytesPerRow + i * 4;
            yuv = rgb2yuv((float)*(pixel) / 255.0, (float)*(pixel + 1) / 255.0, (float)*(pixel + 2) / 255.0);
            radiance_index = (j * width) * 4 + i * 4;
            
            r = radiances[radiance_index];
            g = radiances[radiance_index + 1];
            b = radiances[radiance_index + 2];
            lum = radiances[radiance_index + 3];
            
            
            lum = lum * alpha / lw_global;
            lum = (lum * (1.0 + lum / (l_white * l_white))) / (lum + 1.0);
            
            
            yuv.y = lum;
            rgb = yuv2rgb(yuv.y, yuv.u, yuv.v);
            
            *(pixel) = (int)MAX(MIN(round(rgb.r * 255.0), 255.0), 0.0);
            *(pixel + 1) = (int)MAX(MIN(round(rgb.g * 255.0), 255.0), 0.0);
            *(pixel + 2) = (int)MAX(MIN(round(rgb.b * 255.0), 255.0), 0.0);
            
        }
    }


    int* labels = (int*)malloc(sizeof(int) * width * height);
    float* strength = (float*)malloc(sizeof(float) * width * height);
    int updated = 1;
    int strength_index;

    
    for (j = 0 ; j < height; j++)
    {
        for (i = 0; i < width; i++)
        {
            // ピクセルのポインタを取得する
            pixel = buffer + j * bytesPerRow + i * 4;
            radiance_index = (j * width) * 4 + i * 4;
            
            lum = radiances[radiance_index + 3];
            
            if(lum > lw_global){
                strength[j * width + i] = 1.0;
                labels[j * width + i] = 1;
            }else{
                strength[j * width + i] = 0.0;
                labels[j * width + i] = 0;

            }
        }
    }
    
    int neigh_x, neigh_y;
    
    while (updated > 0) {
        updated = 0;
        for (j = 0 ; j < height; j++)
        {
            for (i = 0; i < width; i++)
            {
                // ピクセルのポインタを取得する
                pixel = buffer + j * bytesPerRow + i * 4;
                radiance_index = (j * width) * 4 + i * 4;
                
                lum = radiances[radiance_index + 3];
                
                for(neigh_y = -1; neigh_y < 2;neigh_y++){
                    if(j + neigh_y < 0 || j + neigh_y >= height){
                        continue;
                    }
                    for(neigh_x = -1; neigh_x < 2;neigh_x++){
                        if(i + neigh_x < 0 || i + neigh_x >= width){
                            continue;
                        }
                        radiance_index = ((j + neigh_y) * width) * 4 + (i + neigh_x) * 4;
                        _lum = radiances[radiance_index + 3];
                        g = 1.0 - _lum / l_white;
                        if(g * strength[(j + neigh_y) * width + (i + neigh_x)] > strength[j * width + i]){
                            strength[j * width + i] = g * strength[(j + neigh_y) * width + (i + neigh_x)];
                            labels[j * width + i] = labels[(j + neigh_y) * width + (i + neigh_x)];
                            updated++;
                        }
                    }
                }
            }
        }
        NSLog(@"updated: %d", updated);
    }
    
    for (j = 0 ; j < height; j++)
    {
        for (i = 0; i < width; i++)
        {
            // ピクセルのポインタを取得する
            pixel = buffer + j * bytesPerRow + i * 4;
            yuv = rgb2yuv((float)*(pixel) / 255.0, (float)*(pixel + 1) / 255.0, (float)*(pixel + 2) / 255.0);
            radiance_index = (j * width) * 4 + i * 4;
            strength_index = j * width + i;
            
            r = radiances[radiance_index];
            g = radiances[radiance_index + 1];
            b = radiances[radiance_index + 2];
            lum = radiances[radiance_index + 3];
            
            
            lum = lum * (1.0 - (strength[strength_index]) / (1.0 + strength[strength_index])) / lw_global;
            lum = (lum * 1.0 + (lum / (l_white * l_white))) / (lum + 1.0);
            
            
            yuv.y = lum;
            rgb = yuv2rgb(yuv.y, yuv.u, yuv.v);
            
            *(pixel) = (int)MAX(MIN(round(rgb.r * 255.0), 255.0), 0.0);
            *(pixel + 1) = (int)MAX(MIN(round(rgb.g * 255.0), 255.0), 0.0);
            *(pixel + 2) = (int)MAX(MIN(round(rgb.b * 255.0), 255.0), 0.0);
            
        }
    }
    
    
     
     for (j = 0 ; j < height; j++)
     {
     for (i = 0; i < width; i++)
     {
     // ピクセルのポインタを取得する
     pixel = buffer + j * bytesPerRow + i * 4;
     yuv = rgb2yuv((float)*(pixel) / 255.0, (float)*(pixel + 1) / 255.0, (float)*(pixel + 2) / 255.0);
     radiance_index = (j * width) * 4 + i * 4;
     
     r = radiances[radiance_index];
     g = radiances[radiance_index + 1];
     b = radiances[radiance_index + 2];
     lum = radiances[radiance_index + 3];
     
     
     lum = lum * alpha / lw_global;
     lum = (lum * 1.0 + (lum / (l_white * l_white))) / (lum + 1.0);
     
     
     yuv.y = lum;
     rgb = yuv2rgb(yuv.y, yuv.u, yuv.v);
     
     *(pixel) = (int)MAX(MIN(round(rgb.r * 255.0), 255.0), 0.0);
     *(pixel + 1) = (int)MAX(MIN(round(rgb.g * 255.0), 255.0), 0.0);
     *(pixel + 2) = (int)MAX(MIN(round(rgb.b * 255.0), 255.0), 0.0);
     
     }
     }

     
     int range = 3;
     int count;
     int tmp_x;
     int tmp_y;
     
     for (j = 0 ; j < height; j++)
     {
     for (i = 0; i < width; i++)
     {
     // ピクセルのポインタを取得する
     lw = 0.0;
     count = 0;
     for(int x = 0;x < range * 2;x++){
     for(int y = 0;y < range * 2;y++){
     tmp_y = (j - (y - range));
     tmp_x = (i - (x - range));
     if(tmp_x < 0){
     continue;
     }
     if(tmp_x > width){
     continue;
     }
     if(tmp_y < 0){
     continue;
     }
     if(tmp_y > height){
     continue;
     }
     radiance_index = tmp_y * width * 4 + tmp_x * 4;
     lw += radiances[radiance_index + 3];
     count++;
     }
     }
     lw /= (float)(count);
     
     pixel = buffer + j * bytesPerRow + i * 4;
     yuv = rgb2yuv((float)*(pixel) / 255.0, (float)*(pixel + 1) / 255.0, (float)*(pixel + 2) / 255.0);
     radiance_index = (j * width) * 4 + i * 4;
     
     r = radiances[radiance_index];
     g = radiances[radiance_index + 1];
     b = radiances[radiance_index + 2];
     lum = radiances[radiance_index + 3];
     
     
     lum = lum * (1.0 - lw / (1.0 + lw)) / lw_global;
     lum = (lum * 1.0 + (lum / (l_white * l_white))) / (lum + 1.0);
     
     
     yuv.y = lum;
     rgb = yuv2rgb(yuv.y, yuv.u, yuv.v);
     
     *(pixel) = (int)MAX(MIN(round(rgb.r * 255.0), 255.0), 0.0);
     *(pixel + 1) = (int)MAX(MIN(round(rgb.g * 255.0), 255.0), 0.0);
     *(pixel + 2) = (int)MAX(MIN(round(rgb.b * 255.0), 255.0), 0.0);
     
     }
     }

    
    for (j = 0 ; j < height; j++)
    {
        for (i = 0; i < width; i++)
        {
            // ピクセルのポインタを取得する
            pixel = buffer + j * bytesPerRow + i * 4;
            yuv = rgb2yuv((float)*(pixel) / 255.0, (float)*(pixel + 1) / 255.0, (float)*(pixel + 2) / 255.0);
            radiance_index = (j * width) * 4 + i * 4;
     
            r = radiances[radiance_index];
            g = radiances[radiance_index + 1];
            b = radiances[radiance_index + 2];
            lum = radiances[radiance_index + 3];
        
            
            lum = lum * alpha / lw;
            lum = (lum * 1.0 + (lum / (l_white * l_white))) / (lum + 1.0);
            

            yuv.y = lum;
            rgb = yuv2rgb(yuv.y, yuv.u, yuv.v);
            
            *(pixel) = (int)MAX(MIN(round(rgb.r * 255.0), 255.0), 0.0);
            *(pixel + 1) = (int)MAX(MIN(round(rgb.g * 255.0), 255.0), 0.0);
            *(pixel + 2) = (int)MAX(MIN(round(rgb.b * 255.0), 255.0), 0.0);
            
        }
    }
    
        for (j = 0 ; j < height; j++)
    {
        for (i = 0; i < width; i++)
        {
            // ピクセルのポインタを取得する
            pixel = buffer + j * bytesPerRow + i * 4;
            radiance_index = (j * width) * 4 + i * 4;
            
            r = radiances[radiance_index];
            g = radiances[radiance_index + 1];
            b = radiances[radiance_index + 2];
            
            r = r * ratio;
            r = MAX(MIN(round(r * 255.0), 255.0), 0.0);
            *(pixel + 0) = (int)r;
            g = g * ratio;
            g = MAX(MIN(round(g * 255.0), 255.0), 0.0);
            *(pixel + 1) = (int)g;
            b = b * ratio;
            b = MAX(MIN(round(b * 255.0), 255.0), 0.0);
            *(pixel + 2) = (int)b;
        
        }
    }
     
     */
    
    

    
    
    free(radiances);
    free(lum_delta);
    
    // do stuff...
    NSTimeInterval timeInterval = [start timeIntervalSinceNow];
    NSLog(@"time:%lf", -timeInterval);

    
    // 効果を与えたデータを作成する
    CFDataRef effectedData;
    effectedData = CFDataCreate(NULL, buffer, CFDataGetLength(mutableData));
    CFRelease(mutableData);
    CFRelease(data);
    
    // 効果を与えたデータプロバイダを作成する
    CGDataProviderRef effectedDataProvider;
    effectedDataProvider = CGDataProviderCreateWithCFData(effectedData);
    
    // 画像を作成する
    CGImageRef effectedCgImage = CGImageCreate(
                                               width, height,
                                               bitsPerComponent, bitsPerPixel, bytesPerRow,
                                               colorSpace, bitmapInfo, effectedDataProvider,
                                               NULL, shouldInterpolate, intent);
    
    
    resultImage = [[UIImage alloc] initWithCGImage:effectedCgImage];
    
    // 作成したデータを解放する
    CGImageRelease(effectedCgImage);
    CFRelease(effectedDataProvider);
    CFRelease(effectedData);
    return resultImage;
}

- (void)didReceiveMemoryWarning
{
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}

@end
