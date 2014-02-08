//
//  ViewController.h
//  HDRtest
//
//  Created by SSC on 2014/01/18.
//  Copyright (c) 2014年 SSC. All rights reserved.
//

#import <UIKit/UIKit.h>
#import <AssetsLibrary/AssetsLibrary.h>
#import "GPUImage.h"
#import "GPUGaussSeidelFilter.h"
#import "GPUFattalFilter.h"
#import "GPUImagePoissonBlendFilter.h"

@interface ViewController : UIViewController <UINavigationControllerDelegate, UIImagePickerControllerDelegate>

@property (nonatomic, strong) UIImage* loadedImage;

@end
