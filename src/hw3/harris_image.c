#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

//int n : jumlah elemen dalam array
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Buat deskriptor fitur untuk indeks dalam gambar
// image im: sumber image
// int i: indeks dalam gambar untuk piksel yang ingin digambarkan
// returns: deskriptor untuk indeks tersebut
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // mengurangi nilai sentral dari neighbors, 
    //untuk mengkompensasi beberapa perubahan eksposur/pencahayaan
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// menandai tempat titik dalam gambar
// image im: gambar untuk menandai
// point p: tempat untuk menandai pada gambar
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// menandai sudut dengan array deskriptor
// image im: gambar untuk menandai
// descriptor *d: sudut pada gambar
// int n: jumlah deskriptor untuk menandai
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// membuat filter Gaussian 1d
// float sigma: standar deviasi Gaussian
// returns: gambar baris tunggal dari filter
image make_1d_gaussian(float sigma)
{
    // TODO: opsional, buat 1d Gaussian yang dapat dipisahkan
    return make_image(1,1,1);
}

// menghaluskan gambar menggunakan filter Gaussian yang dapat dipisahkan
// image im: gambar untuk menghaluskan
// float sigma: std dev untuk Gaussian
// returns: gambar halus
image smooth_image(image im, float sigma)
{
    if(1){
        image g = make_gaussian_filter(sigma);
        image s = convolve_image(im, g, 1);
        free_image(g);
        return s;
    } else {
        // TODO: opsional, gunakan dua konvolusi dengan filter 1d gaussian
        return copy_image(im);
    }
}

// menghitung struktur matriks dari sebuah gambar
// image im: gambar masukan
// float sigma: std dev digunakan untuk menghitung weighted sum
// returns:struktur matriks channel 1: lx^2, channel 2: ly^2, channel 3: lxly
image structure_matrix(image im, float sigma)
{
    // TODO: menghitung struktur matriks untuk im
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image Ix_img = convolve_image(im, gx_filter, 0);
    image Iy_img = convolve_image(im, gy_filter, 0);
    image gs_filter = make_gaussian_filter(sigma);
    image measure_img = make_image(im.w, im.h, 3);
    for (int y=0; y<im.h; y++) {
        for (int x=0; x<im.w; x++) {
            float Ix2 = pow(get_pixel(Ix_img, x, y, 0), 2);
            float Iy2 = pow(get_pixel(Iy_img, x, y, 0), 2);
            float IxIy = get_pixel(Ix_img, x, y, 0) * get_pixel(Iy_img, x, y, 0);
            set_pixel(measure_img, x, y, 0, Ix2);
            set_pixel(measure_img, x, y, 1, Iy2);
            set_pixel(measure_img, x, y, 2, IxIy);
        }
    }
    image S = convolve_image(measure_img, gs_filter, 1);
    free_image(measure_img);
    free_image(gx_filter);
    free_image(gy_filter);
    free_image(Ix_img);
    free_image(Iy_img);
    free_image(gs_filter);
    return S;
}

//  perkirakan sudut sudut setiap piksel dengan struktur matriks S
// image S: struktur matriks untuk gambar
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    //TODO: isi R, "cornerness" untuk setiap piksel menggunakan matriks 
    // gunakan formulasi det(S) - alpha * trace(S)^2, alpha = .06
    for(int y=0; y<S.h; y++) {
        for (int x = 0; x<S.w; x++) {
            // matrix Smatrix = make_matrix(2,2);
            float a = get_pixel(S, x, y, 0);
            float b = get_pixel(S, x, y, 2);
            float c = get_pixel(S, x, y, 2);
            float d = get_pixel(S, x ,y, 1);
            float det = a*d - b*c;
            float trace = a + d;
            float val = det - 0.06*pow(trace, 2);
            set_pixel(R, x, y, 0, val);
        }
    }
    return R;
}

image nms_image(image im, int w)
{
    image r = copy_image(im);
    // TODO: lakukan NMS pada response map
    // untuk setiap piksel dalam gambar:
    //     untuk neighbors dalam w:
    //         jika neighbor response > pixel response:
    //             atur response menjadi sangat rendah (misal: -999999)
    for(int y = 0; y<im.h; y++) {
        for (int x=0; x<im.w; x++) {
            for (int x1=x-w; x1<=x+w; x1++){
                for (int y1=y-w; y1<=y+w; y1++) {
                    if(get_pixel(im, x, y, 0)<get_pixel(im, x1, y1, 0))
                        set_pixel(r, x, y, 0, -999999);
                }
            }
        }
    }
    return r;
}

// Lakukan deteksi sudut harris dan ekstrak fitur dari sudut
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // menghitung struktur matrix
    image S = structure_matrix(im, sigma);

    // memperkirakan sudut
    image R = cornerness_response(S);

    // jalankan NMS pada responses
    image Rnms = nms_image(R, nms);


    //TODO: menghitung jumlah responses terhadap threshold
    int count = 0; 
    for (int y=0; y<Rnms.h; y++) {
        for (int x=0; x<Rnms.w; x++) {
            if (get_pixel(Rnms, x, y, 0)>thresh)
                count+=1;
        }
    }
    int point[count];
    int index = 0;
    for (int y=0; y<Rnms.h; y++) {
        for (int x=0; x<Rnms.w; x++) {
            if (get_pixel(Rnms, x, y, 0)>thresh) {
                point[index] = y*im.w+x;
                index+=1;
            }
        }
    }
    // atur *n sama dengan jumlah sudut pada gambar
    *n = count;
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: isi array *d dengan deskriptor sudut, gunakan describe_index
    for (int i=0; i<count; i++)
        d[i] = describe_index(im, point[i]);
    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// temukan dan gambar sudut pada gambar
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
