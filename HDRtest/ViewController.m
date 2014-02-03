//
//  ViewController.m
//  HDRtest
//
//  Created by SSC on 2014/01/18.
//  Copyright (c) 2014年 SSC. All rights reserved.
//

#import "ViewController.h"

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
        result = [self process:self.loadedImage];

        
    }    else{
        
        GPUDefaultFilter* adjustment = [[GPUDefaultFilter alloc] init];
        
        GPUImagePicture* picture = [[GPUImagePicture alloc] initWithImage:self.loadedImage];
        [picture addTarget:adjustment];
        
        [picture processImage];
        result = [adjustment imageFromCurrentlyProcessedOutput];

    }
    
    UIImageView* view = [[UIImageView alloc] initWithImage:self.loadedImage];
    CGFloat height = self.loadedImage.size.height * 320.0f / self.loadedImage.size.width;
    [view setFrame:CGRectMake(0.0f, 44.0f, 320.0f, height)];
    [self.view addSubview:view];
    
    NSLog(@"%f, %f", height, self.loadedImage.size.width);
    
    CGFloat y = height + 10.0f;
    
    view = [[UIImageView alloc] initWithImage:result];
    height = self.loadedImage.size.height * 320.0f / self.loadedImage.size.width;
    [view setFrame:CGRectMake(0.0f, 44.0f + y, 320.0f, height)];
    [self.view addSubview:view];
}

double rgb2luminance(double red, double green, double blue)
{
    return 0.299 * red + 0.587 * green + 0.114 * blue;
}

double pixel2luminance(unsigned char* pixel)
{
    return rgb2luminance(*(pixel), *(pixel + 1), *(pixel + 2));
}



void multiply(double* r, double* g, double* b)
{
    *r = *r * *r;
    *g = *g * *g;
    *b = *b * *b;
}

void screen2(double* r, double* g, double* b)
{

    *r = 1.0 - (1.0 - *r) * (1.0 - *r);
    *g = 1.0 - (1.0 - *g) * (1.0 - *g);
    *b = 1.0 - (1.0 - *b) * (1.0 - *b);
}

double _log(double v)
{
    return log(v);
    double t = v - 1.0;
    return t - t * t * 0.5 + t * t * t * 0.33333333333333 - t * t * t * t * 0.25;
}

void screen(double* value)
{
    
    double weight = _log(*value);
    weight = weight * weight * 0.05;
    weight = MIN(weight, 1.0);
    weight = 1.0 - weight;
    *value = 1.0 - (1.0 - *value * weight) * (1.0 - *value * weight);
}

void screen_rgb(double* r, double* g, double* b)
{
    double weight = _log(*r);
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
    double r, g, b;
    r = (double)*(rawpixel + 0) * 0.00392156862745;
    g = (double)*(rawpixel + 1) * 0.00392156862745;
    b = (double)*(rawpixel + 2) * 0.00392156862745;
    
    screen_rgb(&r, &g, &b);
    *(rawpixel + 0) = MAX(MIN((int)(r * 255.0), 255), 0);
    *(rawpixel + 1) = MAX(MIN((int)(g * 255.0), 255), 0);
    *(rawpixel + 2) = MAX(MIN((int)(b * 255.0), 255), 0);
}

typedef struct
{
    double r;
    double g;
    double b;
} RGB;

typedef struct
{
    double y;
    double u;
    double v;
} YUV;

typedef struct
{
    int x;
    int y;
    double ev0;
    double ev2;
    double ev4;
    double ev_2;
    double ev_4;
} Position;

typedef struct
{
    double luminance;
    double luminance_multiplied;
    double luminance_screened;
    double delta_lum_multiplied;
    double delta_lum_screened;
} Pixel;

YUV rgb2yuv(double r, double g, double b)
{
    YUV color = {0.0, 0.0, 0.0};
    color.y = 0.299 * r + 0.587 * g + 0.114 * b;
    color.u = -0.169 * r - 0.331 * g + 0.500 * b;
    color.v = 0.500 * r - 0.419 * g - 0.081 * b;
    return color;
}

RGB yuv2rgb(double y, double u, double v)
{
    RGB color = {0.0, 0.0, 0.0};
    color.r = y + 1.402 * v;
    color.g = y - 0.344 * u - 0.714 * v;
    color.b = y + 1.772 * u;
    return color;
}

void calc(unsigned char* rawpixel, Pixel* pixel)
{
    double r, g, b, _r, _g, _b;
    r = (double)*(rawpixel + 0) * 0.00392156862745;
    g = (double)*(rawpixel + 1) * 0.00392156862745;
    b = (double)*(rawpixel + 2) * 0.00392156862745;
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

void inverse(double input[][100])
{
    double inv_a[100][100]; //ここに逆行列が入る
    double buf; //一時的なデータを蓄える
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

double absd(double value)
{
    if(value  < 0.0){
        return -value;
    }
    return value;
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
    
    double r, g, b, _r, _g, _b, ev0, ev2, ev4, ev_2, ev_4, exp0, exp2, exp4, exp_2, exp_4;
    
    exp0 = 1.0;
    exp2 = 0.5;
    exp4 = 0.25;
    exp_2 = 2.0;
    exp_4 = 4.0;
    
    double lum, _lum;
    double tmp;
    double *radiances;
    double *radiances_pixel;
    double *lum_delta;
    int radiance_index;
    
    double max_pixel_value = (1.0 + 1.0 / exp2 + 1.0 / exp4 + 1.0 / exp_2 + 1.0 / exp_4) / 3.0;
    double min_pixel_value = 0.0;
    double lw;
    double lw_global;
    double l_white = 0.0;
    double lw_sum = 0.0;
    
    
    radiances = (double*)malloc(sizeof(double) * width * height * 4);
    lum_delta = (double*)malloc(sizeof(double) * width * height);

    for (j = 0 ; j < height; j++)
    {
        for (i = 0; i < width; i++)
        {
            // ピクセルのポインタを取得する
            pixel = buffer + j * bytesPerRow + i * 4;
            radiance_index = (j * width) * 4 + i * 4;
            //lum = pixel2luminance(pixel) / 255;
            
            
            r = (double)*(pixel + 0) / 255.0;
            g = (double)*(pixel + 1) / 255.0;
            b = (double)*(pixel + 2) / 255.0;

                ev2 = r * r;
                ev4 = ev2 * ev2;
                ev_2 = r;
                screen(&ev_2);
                ev_4 = ev_2;
                screen(&ev_4);
                tmp = r + ev2 / exp2 + ev4 / exp4 + ev_2 / exp_2 + ev_4 / exp_4;
                radiances[radiance_index] = tmp / 5.0;
                ev2 = g * g;
                ev4 = ev2 * ev2;
                ev_2 = g;
                screen(&ev_2);
                ev_4 = ev_2;
                screen(&ev_4);
                tmp = g + ev2 / exp2 + ev4 / exp4 + ev_2 / exp_2 + ev_4 / exp_4;
                radiances[radiance_index + 1] = tmp / 5.0;
                ev2 = b * b;
                ev4 = ev2 * ev2;
                ev_2 = b;
                screen(&ev_2);
                ev_4 = ev_2;
                screen(&ev_4);
                tmp = b + ev2 / exp2 + ev4 / exp4 + ev_2 / exp_2 + ev_4 / exp_4;
                radiances[radiance_index + 2] = tmp / 5.0;
            
            
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
            
            
            r = (double)*(pixel + 0) / 255.0;
            g = (double)*(pixel + 1) / 255.0;
            b = (double)*(pixel + 2) / 255.0;
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
    double alpha = 0.18;
    double ratio = 1.0 / l_white;
    
    double l = 1.0;
    int pyramid_count = 0;
    int base_size = (width > height) ? width : height;
    while (base_size > 128) {
        base_size /= 2;
        pyramid_count++;
    }
    double* gausian_pyramids = (double *)malloc(sizeof(double) * width * height * pyramid_count);
    
    /*
    int cf = (int)pow(2.0, pyramid_count - 1) * 2 + 1;
    double* coefficient = (double*)malloc(sizeof(double) * (cf * (cf + 1)) / 2);
    coefficient[0] = 1.0;
    coefficient[1] = 1.0;
    coefficient[2] = 1.0;
    int index = 0;
    double sum;
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
    int range = 4;
    double N;
    double weight, sum;
    int x, y;
    double* g_table = (double*)malloc(sizeof(double) * width * height);
    double sigma;
    
    for(int n = 0;n < pyramid_count;n++){
        sigma = 0.3 * ((double)(range) / 2.0 - 1.0) + 0.8;
        //sigma = range;
        for (y = 0; y < height; y++) {
            for(x = 0;x < width;x++){
                N = 0.0;
                sum = 0.0;
                for(int xx = -range;xx <= range;xx++){
                    if(x + xx < 0 || x + xx >= width){
                        continue;
                    }
                    _lum = radiances[y * width * 4 + (x + xx) * 4 + 3];
                    weight = exp(-1 * ((xx * xx) / (2.0 * sigma * sigma)));
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
                for(int yy = -range;yy <= range;yy++){
                    if(y + yy < 0 || y + yy >= height){
                        continue;
                    }
                    _lum = g_table[(y + yy) * width + x];
                    weight = exp(-1 * ((yy * yy)) / (2.0 * sigma * sigma));
                    N += weight;
                    sum += _lum * weight;
                }
                gausian_pyramids[width * height * n + y * width + x] = sum / N;
            }
        }
        range *= 2;
        NSLog(@"pyramid");
    }
    
    free(g_table);
    
    
    
    YUV yuv;
    RGB rgb;
    

    double* h_nabla = (double*)malloc(sizeof(double) * width * height * 2);
    double* g_div = (double*)malloc(sizeof(double) * width * height);
    double* phis = (double*)malloc(sizeof(double) * width * height);
    double* result = (double*)malloc(sizeof(double) * width * height);
    double h_x;
    double h_y;
    double h;
    double phi;
    double phi_max = 0.0;
    double L = 0.8;
    
    double g_div_avg = 0.0;
    
    
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
    double beta = 0.8;

    
    for (int j = 1 ; j < height - 1; j++)
    {
        for (int i = 1; i < width - 1; i++)
        {
            phi = 1.0;
            for(int n = pyramid_count - 1;n >= 0;n--){
                h_x = (gausian_pyramids[width * height * n + j * width + (i + 1)] - gausian_pyramids[width * height * n + j * width + (i - 1)]) / pow(2.0, (double)(n + 2));
                h_y = (gausian_pyramids[width * height * n + (j + 1) * width + i] - gausian_pyramids[width * height * n + (j - 1) * width + i]) / pow(2.0, (double)(n + 2));
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
            phis[j * width + i] = phi;
        }
    }
        
    bool* flags = (bool*)malloc(sizeof(bool) * width * height);
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
            //h_nabla[j * width + i + 0] *= phi;
            //h_nabla[j * width + i + 1] *= phi;
            
            g_div[j * width + i] += phi * (h_nabla[j * width * 2 + (i + 1) * 2 + 0] + h_nabla[j * width * 2 + (i - 1) * 2 + 0]);
            g_div[j * width + i] += phi * (h_nabla[(j + 1) * width * 2 + i * 2 + 1] + h_nabla[(j - 1) * width * 2 + i * 2 + 1]);
            
            /*
             g_div[j * width + i] += h_nabla[j * width * 2 + i * 2 + 0] - h_nabla[j * width * 2 + (i - 1) * 2 + 0];
             g_div[j * width + i] += h_nabla[j * width * 2 + i * 2 + 1] - h_nabla[(j - 1) * width * 2 + i * 2 + 1];
             */
        }
    }
    
    double prev_i;
    double threshold = 0.1;
    double max_i = 0.001;
    double current_err = 0.0;
    double max_err = 0.0;
    int updated = 1;
    bool cflag = true;
    
    while (threshold > 0.001) {
        threshold /= 10.0;
        NSLog(@"threshold:%lf", threshold);
        updated = 1;
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
                    if(j == 513 && i == 22){
                        ;
                    }

                    prev_i = result[j * width + i];
                    result[j * width + i] = 0.25 * (result[j * width + i + 1] + result[j * width + i - 1] + result[(j + 1) * width + i] + result[(j - 1) * width + i] - g_div[j * width + i]);
                    result[j * width + i] = prev_i + 1.25 * (result[j * width + i] - prev_i);
                    if(absd(result[j * width + i]) > max_i){
                        max_i = absd(result[j * width + i]);
                    }
                    current_err = absd(prev_i - result[j * width + i]) / max_i;
                    if(current_err > threshold){
                        updated++;
                        flags[j * width + i] = !cflag;
                        flags[j * width + i - 1] = !cflag;
                        flags[j * width + i + 1] = !cflag;
                        flags[(j + 1) * width + i] = !cflag;
                        flags[(j - 1) * width + i] = !cflag;
                    }else{
                        flags[j * width + i] = cflag;
                    }
                }
            }
            NSLog(@"Updated:%d", updated);
        }
    }

    
    double* new_h_nabla = (double*)malloc(sizeof(double) * width * height * 2);
    double* new_g_div = (double*)malloc(sizeof(double) * width * height);
    
    for (int j = 0 ; j < height; j++)
    {
        for (int i = 0; i < width; i++)
        {

            if(j == 2447  && i == 3263){
                ;
            }

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
            
            new_h_nabla[j * width * 2 + i * 2 + 0] = result[j * width + i + 1] - lum;    //nabla Hx
            new_h_nabla[j * width * 2 + i * 2 + 1] = result[(j + 1) * width + i] - lum;    //nabla Hy
            
            
            new_g_div[j * width + i] += new_h_nabla[j * width * 2 + i * 2 + 0] - new_h_nabla[j * width * 2 + (i - 1) * 2 + 0];
            new_g_div[j * width + i] += new_h_nabla[j * width * 2 + i * 2 + 1] - new_h_nabla[(j - 1) * width * 2 + i * 2 + 1];
            
            
            if(i == 22 && j == 513){
                NSLog(@"original g_div:%lf", g_div[j * width + i]);
                NSLog(@"new g_div:%lf", new_g_div[j * width + i]);
                
            }
            
        }
    }
    l_white = 0.0;
    
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

        }
    }
    
    NSLog(@"l_white:%lf", l_white);
    
    ratio = 1.0 / l_white;
    double s = 0.4;
    
    for (int j = 1 ; j < height - 1; j++)
    {
        for (int i = 1; i < width - 1; i++)
        {
            radiance_index = (j * width) * 4 + i * 4;
            pixel = buffer + j * bytesPerRow + i * 4;
            yuv = rgb2yuv((double)*(pixel) / 255.0, (double)*(pixel + 1) / 255.0, (double)*(pixel + 2) / 255.0);
            
            yuv.y = result[j * width + i];
            rgb = yuv2rgb(yuv.y, yuv.u, yuv.v);
            
            if(j == 513 && i == 22){
                NSLog(@"original g_div:%lf", g_div[j * width + i]);
                NSLog(@"new g_div:%lf", new_g_div[j * width + i]);
            }

            tmp = radiances[radiance_index + 0];
            tmp = radiances[radiance_index + 3];
            tmp = result[j * width + i];
                        
            
            rgb.r = pow(radiances[radiance_index + 0] / radiances[radiance_index + 3], s) * result[j * width + i];
            rgb.g = pow(radiances[radiance_index + 1] / radiances[radiance_index + 3], s) * result[j * width + i];
            rgb.b = pow(radiances[radiance_index + 2] / radiances[radiance_index + 3], s) * result[j * width + i];
            
            //rgb.r = rgb.g = rgb.b = exp(gausian_pyramids[width * height * 4 + j * width + i]);
            
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
            yuv = rgb2yuv((double)*(pixel) / 255.0, (double)*(pixel + 1) / 255.0, (double)*(pixel + 2) / 255.0);
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
    double* strength = (double*)malloc(sizeof(double) * width * height);
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
            yuv = rgb2yuv((double)*(pixel) / 255.0, (double)*(pixel + 1) / 255.0, (double)*(pixel + 2) / 255.0);
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
     yuv = rgb2yuv((double)*(pixel) / 255.0, (double)*(pixel + 1) / 255.0, (double)*(pixel + 2) / 255.0);
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
     lw /= (double)(count);
     
     pixel = buffer + j * bytesPerRow + i * 4;
     yuv = rgb2yuv((double)*(pixel) / 255.0, (double)*(pixel + 1) / 255.0, (double)*(pixel + 2) / 255.0);
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
            yuv = rgb2yuv((double)*(pixel) / 255.0, (double)*(pixel + 1) / 255.0, (double)*(pixel + 2) / 255.0);
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
