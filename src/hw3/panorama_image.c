#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"

// membandingkan untuk mencocokkan
// const void *a, *b: pointers untuk membandingkan
// returns: hasil bernilai 0 jika sama, 1 jika a > b, dan -1 jika a < b
int match_compare(const void *a, const void *b)
{
    match *ra = (match *)a;
    match *rb = (match *)b;
    if (ra->distance < rb->distance) return -1;
    else if (ra->distance > rb->distance) return  1;
    else return 0;
}

// fungsi bantuan untuk membuat 2d points
// float x, y: koordinat dari point
// returns: point
point make_point(float x, float y)
{
    point p;
    p.x = x; p.y = y;
    return p;
}

// Tempatkan dua gambar berdampingan di kanvas, untuk menggambar piksel yang cocok
// image a, b: gambar untuk ditempatkan
// returns: gambar a dan b berdampingan
image both_images(image a, image b)
{
    image both = make_image(a.w + b.w, a.h > b.h ? a.h : b.h, a.c > b.c ? a.c : b.c);
    int i,j,k;
    for(k = 0; k < a.c; ++k){
        for(j = 0; j < a.h; ++j){
            for(i = 0; i < a.w; ++i){
                set_pixel(both, i, j, k, get_pixel(a, i, j, k));
            }
        }
    }
    for(k = 0; k < b.c; ++k){
        for(j = 0; j < b.h; ++j){
            for(i = 0; i < b.w; ++i){
                set_pixel(both, i+a.w, j, k, get_pixel(b, i, j, k));
            }
        }
    }
    return both;
}

// menarik garis antara piksel yang cocok dalam dua gambar
// image a, b: dua gambar yang memiliki kecocokan
// match *matches: array yang cocok antara a dan b
// int n: jumlah yang cocok
image draw_matches(image a, image b, match *matches, int n, int inliers)
{
    image both = both_images(a, b);
    int i,j;
    for(i = 0; i < n; ++i){
        int bx = matches[i].p.x; 
        int ex = matches[i].q.x; 
        int by = matches[i].p.y;
        int ey = matches[i].q.y;
        for(j = bx; j < ex + a.w; ++j){
            int r = (float)(j-bx)/(ex+a.w - bx)*(ey - by) + by;
            set_pixel(both, j, r, 0, i<inliers?0:1);
            set_pixel(both, j, r, 1, i<inliers?1:0);
            set_pixel(both, j, r, 2, 0);
        }
    }
    return both;
}

// gambar kecocokan dengan inlier berwarna hijau di antara dua gambar
// image a, b: dua gambar yang cocok
// matches *
image draw_inliers(image a, image b, matrix H, match *m, int n, float thresh)
{
    int inliers = model_inliers(H, m, n, thresh);
    image lines = draw_matches(a, b, m, n, inliers);
    return lines;
}

// temukan sudut, cocokkan, dan gambarkan di antara dua gambar
// image a, b: gambar untuk dicocokkan
// float sigma: gaussian untuk detektor sudut harris
// float thresh: threshold untuk sudut
image find_and_draw_matches(image a, image b, float sigma, float thresh, int nms)
{
    int an = 0;
    int bn = 0;
    int mn = 0;
    descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);
   
    match *m = match_descriptors(ad, an, bd, bn, &mn);

    mark_corners(a, ad, an);
    mark_corners(b, bd, bn);
    image lines = draw_matches(a, b, m, mn, 0);

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);
    return lines;
}

// Menghitung jarak L1 antara array floating point
// float *a, *b: array untuk dibandingkan
// int n: jumlah nilai di setiap array
// returns: l1 jarak antar array (jumlah perbedaan absolut)
float l1_distance(float *a, float *b, int n)
{
    // TODO
    float dis = 0;
    for (int i = 0; i<n; i++) 
        dis += fabs(a[i]-b[i]);
    return dis;
}

// menemukan kecocokan terbaik antara deskriptor dua gambar
// descriptor *a, *b: array deskriptor untuk piksel dalam dua gambar
// int an, bn: jumlah deskriptor dalam array a dan b
match *match_descriptors(descriptor *a, int an, descriptor *b, int bn, int *mn)
{
    int i,j;

    *mn = an;
    match *m = calloc(an, sizeof(match));
    for(j = 0; j < an; ++j){
        // TODO: untuk setiap deskriptor di a, temukan kecocokan terbaik di b
        // record ai sebagai indeks di *a dan bi sebagai indeks di *b
        float min_dist = l1_distance(a[j].data, b[0].data, a[j].n);
        int bind = 0; 
        for(i=1; i<bn; ++i){
            float dist = l1_distance(a[j].data, b[i].data, a[j].n);
            if (dist<min_dist){ 
                bind = i;
                min_dist = dist;
            }
        }       
        m[j].ai = j;
        m[j].bi = bind; // <- index di b
        m[j].p = a[j].p;
        m[j].q = b[bind].p;
        m[j].distance = min_dist; // <- jarak L1 terkecil
    }

    int count = 0;
    
    int *seen = calloc(bn, sizeof(int));
    // TODO: kecocokan satu-ke-satu (injective)
    // urutkan kecocokan berdasarkan jarak dengan match_compare dan qsort
    // arahkan ke elemen yang sama di b
    qsort(m, an, sizeof(match), match_compare);
    for(int i=0; i<an; i++) {
        if (seen[m[i].bi]==0) {
            seen[m[i].bi]=1;
            m[count] = m[i];
            count+=1;
        }
    }
    *mn = count;
    free(seen);
    return m;
}

// terapkan transformasi proyektif ke suatu point
point project_point(matrix H, point p)
{
    matrix c = make_matrix(3, 1);
    // TODO:
    // koordinat homogen setara dengan skalar
    // harus dibagi dengan sesuatu
    c.data[0][0] = p.x;
    c.data[1][0] = p.y;
    c.data[2][0] = 1;
    matrix a = matrix_mult_matrix(H, c);
    point q = make_point(0, 0);
    q.x = a.data[0][0]/a.data[2][0];
    q.y = a.data[1][0]/a.data[2][0];
    return q;
}

// hitung jarak L2 antara dua titik
float point_distance(point p, point q)
{
    // TODO:
    float dist = sqrt(pow(p.x-q.x,2)+pow(p.y-q.y,2));
    return dist;
}

// hitung jumlah inlier
// matrix H: homografi antara sistem koordinat
// match *m: untuk menghitung inlier/outlier
// float thresh: threshold untuk inlier
int model_inliers(matrix H, match *m, int n, float thresh)
{
    int i;
    int count = 0;
    // TODO: hitung jumlah kecocokan yang merupakan inlier
    // yaitu jarak(H*p, q) < thresh
    for(i=0; i<n; i++) {
        point proj_p = project_point(H, m[i].p);
        float dist = point_distance(proj_p, m[i].q);
        
        if (dist<thresh) {
            match temp = m[count];
            m[count]=m[i];
            m[i]=temp;
            count+=1;
        }
    }
    return count;
}

// fungsi random shuffle untuk RANSAC
void randomize_matches(match *m, int n)
{
    // TODO: mengimplementasikan Fisher-Yates untuk shuffle array 
    for(int i=n-1; i>0; i--){
        int j = rand()%(i+1);
        match temp = m[i];
        m[i] = m[j];
        m[j] = temp;
    }
}

// menghitung homografi antara dua gambar yang diberi piksel yang cocok
// match *matches: titik pencocokan antar gambar
// int n: jumlah kecocokan yang digunakan dalam menghitung homografi
// returns: matriks yang merepresentasikan homografi H yang memetakan citra a ke citra b
matrix compute_homography(match *matches, int n)
{
    matrix M = make_matrix(n*2, 8);
    matrix b = make_matrix(n*2, 1);

    int i;
    for(i = 0; i < n; ++i){
        double x  = matches[i].p.x;
        double xp = matches[i].q.x;
        double y  = matches[i].p.y;
        double yp = matches[i].q.y;
        // TODO: isi matriks M dan b
        M.data[2*i][0] = x;
        M.data[2*i][1] = y;
        M.data[2*i][2] = 1;
        M.data[2*i][6] = -x*xp;
        M.data[2*i][7] = -y*xp;
        M.data[2*i+1][3] = x;
        M.data[2*i+1][4] = y;
        M.data[2*i+1][5] = 1;
        M.data[2*i+1][6] = -x*yp;
        M.data[2*i+1][7] = -y*yp;
        b.data[2*i][0] = xp;
        b.data[2*i+1][0] = yp;
    }
    matrix a = solve_system(M, b);
    
    free_matrix(M); free_matrix(b); 

    // Jika solusi tidak ditemukan, kembalikan matriks kosong
    matrix none = {0};
    if(!a.data) return none;
    matrix H = make_matrix(3, 3);
    // TODO: isi homografi H berdasarkan hasil di a
    for(int i=0; i<8; i++) 
        H.data[i/3][i%3] = a.data[i][0];
    H.data[2][2] = 1;
    free_matrix(a);
    return H;
}

// lakukan konsensus sampel acak untuk menghitung homografi
// int k: jumlah iterasi yang akan dijalankan
// int cutoff: cutoff inlier untuk keluar lebih awal
matrix RANSAC(match *m, int n, float thresh, int k, int cutoff)
{
    int best = 0;
    matrix Hb = make_translation_homography(256, 0);
    // TODO: algoritma RANSAC 
    // for k iterations:
    //     shuffle nilai kecocokan
    //     menghitung homografi dengan beberapa kecocokan 
    //     jika homografi baru lebih baik daripada yang lama:
    //         menghitung homografi yang diperbarui menggunakan semua inlier
    //         jika itu lebih baik dari cutoff:
    //             kembalikan
    for(int i=0; i<k; i++) {
        randomize_matches(m, n);
        matrix Hi = compute_homography(m, 4);
        int inliers = model_inliers(Hi, m, n, thresh);
        if (inliers>best) {
            matrix Hi2 = compute_homography(m, inliers);
            best = inliers;
            Hb = Hi2;
            if (model_inliers(Hb, m, n, thresh)>=cutoff)
                return Hb;
        }
    }
    return Hb;
}

// menggabungkan dua gambar menggunakan transformasi proyektif
// image a, b: gambar untuk digabungkan
// matrix H: homografi dari koordinat citra a ke koordinat citra b
// returns: gabungan dari beberapa gambar 
image combine_images(image a, image b, matrix H)
{
    matrix Hinv = matrix_invert(H);

    // proyeksikan sudut-sudut gambar b ke dalam koordinat gambar a
    point c1 = project_point(Hinv, make_point(0,0));
    point c2 = project_point(Hinv, make_point(b.w-1, 0));
    point c3 = project_point(Hinv, make_point(0, b.h-1));
    point c4 = project_point(Hinv, make_point(b.w-1, b.h-1));

    // temukan sudut kiri atas dan kanan bawah gambar b yang dilengkungkan menjadi gambar a
    point topleft, botright;
    botright.x = MAX(c1.x, MAX(c2.x, MAX(c3.x, c4.x)));
    botright.y = MAX(c1.y, MAX(c2.y, MAX(c3.y, c4.y)));
    topleft.x = MIN(c1.x, MIN(c2.x, MIN(c3.x, c4.x)));
    topleft.y = MIN(c1.y, MIN(c2.y, MIN(c3.y, c4.y)));

    // temukan seberapa besar seharusnya gambar baru dan offset gambar a
    int dx = MIN(0, topleft.x);
    int dy = MIN(0, topleft.y);
    int w = MAX(a.w, botright.x) - dx;
    int h = MAX(a.h, botright.y) - dy;

    if(w > 7000 || h > 7000){
        fprintf(stderr, "output too big, stopping\n");
        return copy_image(a);
    }

    int i,j,k;
    image c = make_image(w, h, a.c);
    
    // paste gambar a ke dalam gambar baru yang diimbangi dengan dx dan dy
    for(k = 0; k < a.c; ++k){
        for(j = 0; j < a.h; ++j){
            for(i = 0; i < a.w; ++i){
                // TODO: 
                float val = get_pixel(a, i, j, k);
                set_pixel(c, i-dx, j-dy, k, val);
            }
        }
    }

    // TODO: paste di gambar b juga  
    for(int k=0; k<a.c; k++) {
        for (int y=topleft.y; y<botright.y; y++) {
            for (int x=topleft.x; x<botright.x; x++) {
                point proj = project_point(H, make_point(x, y));
                if(proj.x<b.w && proj.x>=0 && proj.y<b.h && proj.y>=0){
                    float val = bilinear_interpolate(b, proj.x, proj.y, k);
                    set_pixel(c, x-dx, y-dy, k, val);
                }
            }
        }
    }
    
    return c;
}

// buat panorama antara dua gambar
// image a, b: gambar untuk digabungkan
// float sigma: gaussian untuk detektor sudut harris
// float thresh: threshold untuk sudut/tidak ada sudut
// float inlier_thresh: threshold untuk RANSAC inliers
image panorama_image(image a, image b, float sigma, float thresh, int nms, float inlier_thresh, int iters, int cutoff)
{
    srand(10);
    int an = 0;
    int bn = 0;
    int mn = 0;
    
    // hitung sudut dan deskriptor
    descriptor *ad = harris_corner_detector(a, sigma, thresh, nms, &an);
    descriptor *bd = harris_corner_detector(b, sigma, thresh, nms, &bn);

    // cari kecocokan
    match *m = match_descriptors(ad, an, bd, bn, &mn);

    // jalankan RANSAC untuk mendapatkan homography
    matrix H = RANSAC(m, mn, inlier_thresh, iters, cutoff);

    if(1){
        // tandai sudut dan cocokkan antara gambar
        mark_corners(a, ad, an);
        mark_corners(b, bd, bn);
        image inlier_matches = draw_inliers(a, b, H, m, mn, inlier_thresh);
        save_image(inlier_matches, "inliers");
    }

    free_descriptors(ad, an);
    free_descriptors(bd, bn);
    free(m);

    // gabungkan gambar dengan homografi
    image comb = combine_images(a, b, H);
    return comb;
}

image cylindrical_project(image im, float f)
{
    //TODO: 
    // menentukan ukuran gambar silinder
    int xc = im.w/2;
    int yc = im.h/2;
    int flat_w = 2*f*atan2(xc/f, 1);
    int flat_h = im.h;
    
    image flat = make_image(flat_w, flat_h, im.c);
    // untuk setiap piksel, lakukan interpolasi
    for (int k=0; k<im.c; k++) {
        for (int y=0; y<flat.h; y++){
            for (int x=0; x<flat.w; x++) {
                float x2=tan((x-flat.w/2)/f)*f+xc;
                float y2=(y-flat.h/2)/cos((x-flat.w/2)/f)+yc;
                float val = 0;
                if (x2>=0 && x2<im.w && y2>=0 && y2<=im.h)
                    val = bilinear_interpolate(im, x2, y2, k);
                set_pixel(flat, x, y, k, val); 
            }
        }
    }
    return flat;
}
image spherical_project(image im, float f)
{
    //TODO: 
    // menentukan ukuran gambar silinder
    int xc = im.w/2;
    int yc = im.h/2;
    int flat_w = 2*f*atan2(xc/f, 1);
    int flat_h = 2*f*atan2(yc/f, 1);
    
    image flat = make_image(flat_w, flat_h, im.c);
    // untuk setiap piksel, lakukan interpolasi
    for (int k=0; k<im.c; k++) {
        for (int y=0; y<flat.h; y++){
            for (int x=0; x<flat.w; x++) {
                float x2=f*tan((x-flat.w/2)/f)+xc;
                float y2=f*tan((y-flat.h/2)/f)/cos((x-flat.w/2)/f)+yc;
                float val = 0;
                if (x2>=0 && x2<im.w && y2>=0 && y2<=im.h)
                    val = bilinear_interpolate(im, x2, y2, k);
                set_pixel(flat, x, y, k, val); 
            }
        }
    }
    return flat;
}
