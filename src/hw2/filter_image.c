#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

void l1_normalize(image im)
{
    //normalisasi gambar untuk digabungkan menjadi 1
    float sum = 0;
    for (int c = 0; c<im.c; c++) {
        for(int y = 0; y<im.h; y++) {
            for(int x=0; x<im.w; x++) {
                sum+=get_pixel(im, x,y,c);
            }
        }
    }
    for (int c = 0; c<im.c; c++) {
        for(int y = 0; y<im.h; y++) {
            for(int x=0; x<im.w; x++) {
                set_pixel(im, x,y,c, get_pixel(im, x,y,c)/sum);
            }
        }
    }

}

image make_box_filter(int w)
{
    //buat filter kotak persegi dengan ukuran w dengan satu channel 
    image filter = make_image(w,w,1);
    for(int y = 0; y<w; y++) {
        for(int x=0; x<w; x++) {
            set_pixel(filter, x,y,0, 1);
        }
    }
    l1_normalize(filter);
    return filter;
}

image convolve_image(image im, image filter, int preserve)
{
    int half_size = (filter.w-1)/2;
    assert(filter.c==1 || im.c ==filter.c);
    image output_img;
    // jika perserve = 1, maka channel keluaran gambar sama dengan im 
    if(preserve==1) {
        output_img = make_image(im.w, im.h, im.c);
        //looping bagian luar untuk mendapatkan setiap nilai piksel keluaran gambar
        for (int c = 0; c<im.c; c++) {
            for(int y = 0; y<im.h; y++) {
                for(int x=0; x<im.w; x++) {
                    //inside loop : untuk setiap nilai pixel keluaran gambar, kita perlu melakukan convolution
                    float pixel_val = 0;
                    for(int y_conv = y-half_size; y_conv<=y+half_size; y_conv++){
                        for(int x_conv = x-half_size; x_conv<=x+half_size; x_conv++) {
                            float val = get_pixel(im, x_conv, y_conv, c);
                            //periksa jumlah channel filter
                            int filter_c = (filter.c ==1)?0:c;
                            float filter_val = get_pixel(filter, x_conv-x+half_size, y_conv-y+half_size, filter_c);
                            pixel_val += val*filter_val;                            
                        }
                    }
                    set_pixel(output_img, x,y,c,pixel_val);
                }
            }
        }
    }
    //jika perserve = 0, kita perlu menjumlahkan dimensi spasial dan channel untuk menghasilkan channel gambar 1 
    else if (preserve==0) {
        output_img = make_image(im.w, im.h, 1);
        for (int y = 0; y<im.h; y++) {
            for(int x = 0; x<im.w; x++) {
                float sumVal = 0;
                for(int c=0; c<im.c; c++) {
                    float pixel_val = 0;
                    for(int y_conv = y-half_size; y_conv<=y+half_size; y_conv++){
                        for(int x_conv = x-half_size; x_conv<=x+half_size; x_conv++) {
                            float val =get_pixel(im, x_conv,y_conv,c);
                            int filter_c = (filter.c == 1)?0:c;
                            pixel_val += val*get_pixel(filter, x_conv-x+half_size, y_conv-y+half_size, filter_c);
                        }
                    }
                    sumVal += pixel_val;
                }
                set_pixel(output_img, x,y,0,sumVal);
            }
        }
    }
    return output_img;
}

image make_highpass_filter()
{
    /* Bentuk filter : 
0  -1  0
     -1   4 -1
      0  -1  0
    */
    image hp = make_image(3,3,1);
    hp.data[0] = 0;
    hp.data[1] = -1;
    hp.data[2] = 0;
    hp.data[3] = -1;
    hp.data[4] = 4;
    hp.data[5] = -1;
    hp.data[6] = 0;
    hp.data[7] = -1;
    hp.data[8] = 0;
    return hp;
}

image make_sharpen_filter()
{
    /* Bentuk filter :
       0  -1  0
      -1   5 -1
       0  -1  0
    */

    image sp = make_image(3,3,1);
    sp.data[0] = 0;
    sp.data[1] = -1;
    sp.data[2] = 0;
    sp.data[3] = -1;
    sp.data[4] = 5;
    sp.data[5] = -1;
    sp.data[6] = 0;
    sp.data[7] = -1;
    sp.data[8] = 0;
    return sp;
}

image make_emboss_filter()
{
    /* Bentuk filter :
       -2 -1  0
       -1  1  1
        0  1  2
    */
    image eb = make_image(3,3,1);
    eb.data[0] = -2;
    eb.data[1] = -1;
    eb.data[2] = 0;
    eb.data[3] = -1;
    eb.data[4] = 1;
    eb.data[5] = 1;
    eb.data[6] = 0;
    eb.data[7] = 1;
    eb.data[8] = 2;
    return eb;
}

image make_gaussian_filter(float sigma)
{
    // TODO
    int filter_size = (int)((int)ceil(sigma*6)%2==0)?(int)ceil(sigma*6)+1:(int)ceil(sigma*6);
    image gs_filter = make_image(filter_size, filter_size, 1);
    for (int j=0; j<gs_filter.h; j++) {
        for(int i=0; i<gs_filter.w; i++) {
            int x = i - (filter_size-1)/2;
            int y = j - (filter_size-1)/2;
            float val = exp(-(pow(x,2)+pow(y,2))/(2*pow(sigma,2)))/(TWOPI*pow(sigma,2));
            set_pixel(gs_filter, i, j, 0,val);
        }
    }
    l1_normalize(gs_filter);
    return gs_filter;
}

image add_image(image a, image b)
{
    // TODO
    assert(a.w==b.w && a.h==b.h && a.c == b.c);
    image add = make_image(a.w, a.h, a.c);
    for(int c = 0; c<a.c; c++) {
        for(int y = 0; y<a.h; y++) {
            for (int x =0; x<a.w; x++) {
                float val_a = get_pixel(a, x, y, c);
                float val_b = get_pixel(b, x, y, c);
                set_pixel(add, x, y, c, val_a+val_b);
            }
        }
    }
    
    return add;
}


image sub_image(image a, image b)
{
    // TODO
    assert(a.w==b.w && a.h==b.h && a.c == b.c);
    image sub = make_image(a.w, a.h, a.c);
    for(int c = 0; c<a.c; c++) {
        for(int y = 0; y<a.h; y++) {
            for (int x =0; x<a.w; x++) {
                float val_a = get_pixel(a, x, y, c);
                float val_b = get_pixel(b, x, y, c);
                set_pixel(sub, x, y, c, val_a-val_b);
            }
        }
    }
    return sub;
}

image make_gx_filter()
{
    /* Bentuk filter:
       -1  0  1
       -2  0  2
    */ -1  0  1
    image gx = make_image(3,3,1);
    gx.data[0] = -1;
    gx.data[1] = 0;
    gx.data[2] = 1;
    gx.data[3] = -2;
    gx.data[4] = 0;
    gx.data[5] = 2;
    gx.data[6] = -1;
    gx.data[7] = 0;
    gx.data[8] = 1;
    return gx;
}

image make_gy_filter()
{
    /*Bentuk filter:
      -1  -2  -1
       0   0   0
       1   2   1
    */
    image gy = make_image(3,3,1);
    gy.data[0] = -1;
    gy.data[1] = -2;
    gy.data[2] = -1;
    gy.data[3] = 0;
    gy.data[4] = 0;
    gy.data[5] = 0;
    gy.data[6] = 1;
    gy.data[7] = 2;
    gy.data[8] = 1;
    return gy;
}

void feature_normalize(image im)
{
    // TODO
    float min_val = get_pixel(im, 0, 0, 0);
    float max_val = min_val;
    for (int c = 0; c<im.c; c++) {
        for (int y = 0; y<im.h; y++) {
            for (int x = 0; x<im.w; x++) {
                float val = get_pixel(im, x, y, c);
                min_val = (val<min_val)?val:min_val;
                max_val = (val>max_val)?val:max_val;
            }
        }
    }
    float range = max_val - min_val;
    for (int c = 0; c<im.c; c++) {
        for (int y = 0; y<im.h; y++) {
            for (int x = 0; x<im.w; x++) {
                set_pixel(im, x, y, c, (get_pixel(im, x, y, c)-min_val)/range);
            }
        }
    }
    
}

image *sobel_image(image im)
{
    // TODO
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    
    image *r = (image *)calloc(2, sizeof(image));

    image gx = convolve_image(im, gx_filter, 0);
    image gy = convolve_image(im, gy_filter, 0);
    image mag = make_image(im.w, im.h, 1);
    image dir = make_image(im.w, im.h, 1);
    for(int y = 0; y<im.h; y++) {
        for (int x = 0; x<im.w; x++) {
            float mag_val = sqrt(pow(get_pixel(gx, x, y, 0),2) + pow(get_pixel(gy, x, y, 0),2));
            float theta = atan2(get_pixel(gy, x, y, 0),get_pixel(gx, x, y, 0));
            set_pixel(mag, x, y, 0, mag_val);
            set_pixel(dir, x, y, 0, theta);
        }
    }
    r[0] = mag;
    r[1] = dir;
    return r;    
}

image colorize_sobel(image im)
{
    // TODO
    image *res = sobel_image(im);
    image mag = res[0];
    feature_normalize(mag);
    image angle = res[1];
    feature_normalize(angle);
    image im_hsv = make_image(im.w, im.h, 3);
    for (int c = 0; c<3; c++) {
        for (int y = 0; y<im.h; y++) {
            for (int x = 0; x<im.w; x++){
                float h = get_pixel(angle, x, y, 0);
                float s = get_pixel(mag, x, y, 0);
                float v = get_pixel(mag, x, y, 0);
                set_pixel(im_hsv, x, y, 0, h);
                set_pixel(im_hsv, x, y, 1, s);
                set_pixel(im_hsv, x, y, 2, v);
            }
        }
    }
    hsv_to_rgb(im_hsv);
    return im_hsv;
}
