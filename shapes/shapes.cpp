#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

template<class T>
struct vec3 {
    array<T,3> d;
};

template<class T>
vec3<T> operator+(const vec3<T> &a, const vec3<T> &b) {
    vec3<T> ret;
    for (int i = 0; i < 3; ++i) {
        ret.d[i] = a.d[i] + b.d[i];
    }
    return ret;
}

template<class T>
vec3<T> operator-(const vec3<T> &a, const vec3<T> &b) {
    vec3<T> ret;
    for (int i = 0; i < 3; ++i) {
        ret.d[i] = a.d[i] - b.d[i];
    }
    return ret;
}

template<class T>
vec3<T> operator-(const vec3<T> &a) {
    vec3<T> ret;
    for (int i = 0; i < 3; ++i) {
        ret.d[i] = - a.d[i];
    }
    return ret;
}

template<class T>
T operator*(const vec3<T> &a, const vec3<T> &b) {
    T ret = 0;
    for (int i = 0; i < 3; ++i) {
        ret += a.d[i] * b.d[i];
    }
    return ret;
}

template<class T>
vec3<T> operator*(const vec3<T> &a, T b) {
    vec3<T> ret;
    for (int i = 0; i < 3; ++i) {
        ret.d[i] = a.d[i] * b;
    }
    return ret;
}

template<class T>
vec3<T> operator*(T a, const vec3<T> &b) {
    vec3<T> ret;
    for (int i = 0; i < 3; ++i) {
        ret.d[i] = a * b.d[i];
    }
    return ret;
}

template<class T>
bool operator!=(const vec3<T> &a, const vec3<T> &b) {
    return (a.d[0] != b.d[0]) || (a.d[1] != b.d[1]) || (a.d[2] != b.d[2]);
}

template<class T>
vec3<T> cross(const vec3<T> &a, const vec3<T> &b) {
    vec3<T> ret;
    ret.d[0] = a.d[1] * b.d[2] - a.d[2] * b.d[1];
    ret.d[1] = a.d[2] * b.d[0] - a.d[0] * b.d[2];
    ret.d[2] = a.d[0] * b.d[1] - a.d[1] * b.d[0];
    return ret;
}

double abs(const vec3<double> &a) {
    return sqrt(a * a);
}

vec3<double> vec3_int2dbl(const vec3<int> &a) {
    vec3<double> ret;
    ret.d[0] = a.d[0];
    ret.d[1] = a.d[1];
    ret.d[2] = a.d[2];
    return ret;
}

struct plane {
    array<array<array<bool,2>,2>,2> b;
    vec3<double> n;
    int d;
};

struct triangle {
    array<int,3> v;
    vec3<double> n;
    double d;
    array<vec3<double>,3> nn;
    array<double,3> dd;
};

inline bool operator==(const triangle &a, const triangle &b) {
    for (int i = 0; i < 3; ++i) {
        if ((a.v[i] != b.v[i]) || (a.nn[i] != b.nn[i]) || (a.dd[i] != b.dd[i])) {
            return false;
        }
    }
    if ((a.n != b.n) || (a.d != b.d)) {
        return false;
    }
    return true;
}

struct rectangle {
    array<int,4> v;
    vec3<double> n;
    double d;
    array<vec3<double>,4> nn;
    array<double,4> dd;
};

inline bool operator==(const rectangle &a, const rectangle &b) {
    for (int i = 0; i < 4; ++i) {
        if ((a.v[i] != b.v[i]) || (a.nn[i] != b.nn[i]) || (a.dd[i] != b.dd[i])) {
            return false;
        }
    }
    if ((a.n != b.n) || (a.d != b.d)) {
        return false;
    }
    return true;
}

struct shape {
    int bits;
    array<array<array<bool,2>,2>,2> b;
    vector<int> tr;
    vector<int> rc;
    //vector<triangle> tr;
    //vector<rectangle> rc;
};

struct ray {
    vec3<double> o;
    vec3<double> r;
};

// return true if shape has volume > 0, otherwise false
bool check_shape_volume(const shape sh) {
    // number of bits must be at least 4
    int num_bits = 0;
    for (int ix = 0; ix < 2; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int iz = 0; iz < 2; ++iz) {
                if (sh.b[ix][iy][iz]) {
                    ++num_bits;
                }
            }
        }
    }
    if (num_bits < 4) {
        return false;
    }
    
    // check x planes for at least one enabled vertex
    for (int ix = 0; ix < 2; ++ix) {
        bool tmp = false;
        for (int iy = 0; iy < 2; ++iy) {
            for (int iz = 0; iz < 2; ++iz) {
                if (sh.b[ix][iy][iz]) {
                    tmp = true;
                    break;
                }
            }
            if (tmp) {
                break;
            }
        }
        if (!tmp) {
            return false;
        }
    }
    
    // check y planes for at least one enabled vertex
    for (int iy = 0; iy < 2; ++iy) {
        bool tmp = false;
        for (int ix = 0; ix < 2; ++ix) {
            for (int iz = 0; iz < 2; ++iz) {
                if (sh.b[ix][iy][iz]) {
                    tmp = true;
                    break;
                }
            }
            if (tmp) {
                break;
            }
        }
        if (!tmp) {
            return false;
        }
    }
    
    // check z planes for at least one enabled vertex
    for (int iz = 0; iz < 2; ++iz) {
        bool tmp = false;
        for (int ix = 0; ix < 2; ++ix) {
            for (int iy = 0; iy < 2; ++iy) {
                if (sh.b[ix][iy][iz]) {
                    tmp = true;
                    break;
                }
            }
            if (tmp) {
                break;
            }
        }
        if (!tmp) {
            return false;
        }
    }
    
    return true;
}

void get_coordinates(int v, vec3<int> &p) {
    p.d[0] = v & 1;
    p.d[1] = (v & 2) >> 1;
    p.d[2] = (v & 4) >> 2;
}
void get_coordinates(int v, vec3<double> &p) {
    p.d[0] = v & 1;
    p.d[1] = (v & 2) >> 1;
    p.d[2] = (v & 4) >> 2;
}

void print_bit_field(int i, int num_bits) {
    for (int j = num_bits-1; j >= 0; --j) {
        cout << ((i & (1 << j)) >> j);
    }
}

// return distance of r.o to plane, or -1 if no intersection
double ray_plane_intersect(const ray &r, const vec3<double> n, double d) {
    // get ray parameter for plane intersection
    double n_r = n * r.r;
    if (n_r >= 0) {
        // plane parallel to ray or facing in wrong direction
        return -1;
    }
    double t = - (n * r.o + d) / n_r;
    if (t < 0) {
        // plane behind ray origin
        return -1;
    }
    return t;
}

// return distance of r.o to triangle, or -1 if no intersection
double ray_triangle_intersect(const ray &r, const triangle &tr) {
    // ray-plane intersection
    double t = ray_plane_intersect(r, tr.n, tr.d);
    if (t < 0) {
        return -1;
    }
    
    // get intersection point
    vec3<double> p = r.o + t * r.r;
    
    vec3<double> v[3];
    for (int i = 0; i < 3; ++i) {
        get_coordinates(tr.v[i], v[i]);
    }
    
    for (int i = 0; i < 3; ++i) {
        int a = i;
        int b = (i + 1) % 3;
        int c = (i + 2) % 3;

        // edges
        vec3<double> ab = v[b] - v[a];
        vec3<double> ac = v[c] - v[a];
        
        // get plane normal to triangle through edge ab
        vec3<double> n = ac - ((ab * ac) / (ab * ab)) * ab;
        double d = - (n * v[a]);
        
        // no intersection if point is on negative side
        if ((n * p + d) < 0) {
            return -1;
        }
    }
    
    return t;
}

// return distance of r.o to rectangle, or -1 if no intersection
double ray_rectangle_intersect(const ray &r, const rectangle &rc) {
    // ray-plane intersection
    double t = ray_plane_intersect(r, rc.n, rc.d);
    if (t < 0) {
        return -1;
    }
    
    // get intersection point
    vec3<double> p = r.o + t * r.r;
    
    vec3<double> v[4];
    for (int i = 0; i < 4; ++i) {
        get_coordinates(rc.v[i], v[i]);
    }
    
    int ivt[4] = {
        0, 1, 3, 2
    };
    
    for (int i = 0; i < 4; ++i) {
        int a = ivt[i];
        int b = ivt[(i+1)%4];
        int c = ivt[(i+2)%4];
        
        // edges
        vec3<double> ab = v[b] - v[a];
        vec3<double> ac = v[c] - v[a];
        
        // get plane normal to rectangle through edge ab
        vec3<double> n = ac - ((ab * ac) / (ab * ab)) * ab;
        double d = - (n * v[a]);
        
        // no intersection if point is on negative side
        if ((n * p + d) < 0) {
            return -1;
        }
    }

    return t;
}

shape shapes[256];
vector<plane> planes;

int main(int argc, char **argv) {
    // generate all possible shapes and planes
    for (int i = 0; i < 256; ++i) {
        // set vertices
        shapes[i].bits = i;
        for (int j = 0; j < 8; ++j) {
            bool bit = (((i & (1 << j)) >> j) == 1);
            shapes[i].b[j&1][(j&2)>>1][(j&4)>>2] = bit;
        }
        
        /*// set shapes with no volume to empty
        if (!check_shape_volume(shapes[i])) {
            shapes[i].bits = 0;
            for (int j = 0; j < 8; ++j) {
                bool bit = (((i & (1 << j)) >> j) == 1);
                shapes[i].b[j&1][(j&2)>>1][(j&4)>>2] = bit;
            }
        }*/
    }
    
    // get all planes
    plane pl;
    vec3<int> p0, p1, p2, p3;
    for (int v0 = 0; v0 < 6; ++v0) {
        get_coordinates(v0, p0);
        for (int v1 = v0 + 1; v1 < 7; ++v1) {
            get_coordinates(v1, p1);
            for (int v2 = v1 + 1; v2 < 8; ++v2) {
                get_coordinates(v2, p2);
                
                // construct plane equation: nx * x + ny * y + nz * z + d = 0
                vec3<int> n = cross(p2 - p0, p1 - p0);
                int d = - (n * p0);
                
                // check all vertices
                for (int v3 = 0; v3 < 8; ++v3) {
                    get_coordinates(v3, p3);
                    
                    if (n * p3 + d == 0) {
                        pl.b[v3&1][(v3&2)>>1][(v3&4)>>2] = true;
                    } else {
                        pl.b[v3&1][(v3&2)>>1][(v3&4)>>2] = false;
                    }
                }
                
                // add plane if it does not exist already
                bool found = false;
                for (plane pl2 : planes) {
                    bool same = true;
                    for (int i = 0; i < 8; ++i) {
                        if (pl2.b[i&1][(i&2)>>1][(i&4)>>2] != pl.b[i&1][(i&2)>>1][(i&4)>>2]) {
                            same = false;
                            break;
                        }
                    }
                    if (same) {
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    pl.n = vec3_int2dbl(n);
                    pl.d = d;
                    planes.push_back(pl);
                }
            }
        }
    }
    
    // get list of rectangles/triangles for each shape
    std::vector<triangle> triangles;
    std::vector<rectangle> rectangles;
    for (int i = 0; i < 256; ++i) {
        if (shapes[i].bits == 0) {
            continue;
        }
        
        for (int j = 0; j < planes.size(); ++j) {
            // collect points which are enabled in plane and shape
            vector<int> v;
            for (int k = 0; k < 8; ++k) {
                if (shapes[i].b[k&1][(k&2)>>1][(k&4)>>2] && planes[j].b[k&1][(k&2)>>1][(k&4)>>2]) {
                    v.push_back(k);
                }
            }
            
            vec3<double> vv[v.size()];
            for (int k = 0; k < v.size(); ++k) {
                get_coordinates(v[k], vv[k]);
            }
            
            if (v.size() == 3) {
                triangle tr;
                tr.v[0] = v[0];
                tr.v[1] = v[1];
                tr.v[2] = v[2];
                tr.n = planes[j].n;
                tr.d = planes[j].d;

                for (int k = 0; k < 3; ++k) {
                    int a = k;
                    int b = (k + 1) % 3;
                    int c = (k + 2) % 3;
                    
                    // edges
                    vec3<double> ab = vv[b] - vv[a];
                    vec3<double> ac = vv[c] - vv[a];
                    
                    // get plane normal to triangle through edge ab
                    tr.nn[k] = ac * (ab * ab) - (ab * ac) * ab;
                    tr.dd[k] = - (tr.nn[k] * vv[a]);
                }
                
                for (int l = 0; l < 2; ++l) {
                    int k;
                    for (k = 0; k < triangles.size(); ++k) {
                        if (triangles[k] == tr) {
                            break;
                        }
                    }
                    if (k == triangles.size()) {
                        triangles.push_back(tr);
                    }
                    shapes[i].tr.push_back(k);
                
                    // flip orientation
                    tr.n = - tr.n;
                    tr.d = - tr.d;
                }
            } else if (v.size() == 4) {
                rectangle rc;
                rc.v[0] = v[0];
                rc.v[1] = v[1];
                rc.v[2] = v[2];
                rc.v[3] = v[3];
                rc.n = planes[j].n;
                rc.d = planes[j].d;
                
                int ivt[4] = {
                    0, 1, 3, 2
                };
                for (int k = 0; k < 4; ++k) {
                    int a = ivt[k];
                    int b = ivt[(k+1)%4];
                    int c = ivt[(k+2)%4];
                    
                    // edges
                    vec3<double> ab = vv[b] - vv[a];
                    vec3<double> ac = vv[c] - vv[a];
                    
                    // get plane normal to rectangle through edge ab
                    rc.nn[k] = ac * (ab * ab) - (ab * ac) * ab;
                    rc.dd[k] = - (rc.nn[k] * vv[a]);
                }
                
                for (int l = 0; l < 2; ++l) {
                    int k;
                    for (k = 0; k < rectangles.size(); ++k) {
                        if (rectangles[k] == rc) {
                            break;
                        }
                    }
                    if (k == rectangles.size()) {
                        rectangles.push_back(rc);
                    }
                    shapes[i].rc.push_back(k);
                    
                    // flip orientation
                    rc.n = - rc.n;
                    rc.d = - rc.d;
                }
            }
        }
    }
    
    // remove hidden faces by testing ray intersections (4 parallel rays from each direction +/- x, +/- y, +/- z)
    for (int k = 0; k < 256; ++k) {
        bool keep_tr[shapes[k].tr.size()] = { false };
        bool keep_rc[shapes[k].rc.size()] = { false };
        
        for (int dir = 0; dir < 3; ++dir) {
            int od[2] = {
                (dir + 1) % 3, (dir + 2) % 3
            };
        
            for (int pm = -1; pm <= 1; pm += 2) {
                for (int i = 0; i < 2; ++i) {
                    for (int j = 0; j < 2; ++j) {
                        ray r;
                        r.o.d[dir  ] = 1.5 * (pm + 1) - 1;
                        r.o.d[od[0]] = 0.2282 + i * 0.5523;
                        r.o.d[od[1]] = 0.2412 + j * 0.5012;
                        r.r.d[dir  ] = - pm;
                        r.r.d[od[0]] = 0;
                        r.r.d[od[1]] = 0;
                    
                        // check all faces of shape
                        int l0 = -1;
                        double t0 = 1e99;
                        bool tri = false;
                        bool rect = false;
                        for (int l = 0; l < shapes[k].tr.size(); ++l) {
                            double t = ray_triangle_intersect(r, triangles[shapes[k].tr[l]]);
                            if ((t >= 0) && (t < t0)) {
                                tri = true;
                                l0 = l;
                                t0 = t;
                            }
                        }
                        for (int l = 0; l < shapes[k].rc.size(); ++l) {
                            double t = ray_rectangle_intersect(r, rectangles[shapes[k].rc[l]]);
                            if ((t >= 0) && (t <= t0)) {
                                tri  = false;
                                rect = true;
                                l0 = l;
                                t0 = t;
                            }
                        }
                        if (tri) {
                            keep_tr[l0] = true;
                        } else if (rect) {
                            keep_rc[l0] = true;
                        }
                    }
                }
            }
        }
        
        vector<int> tr_tmp = shapes[k].tr;
        vector<int> rc_tmp = shapes[k].rc;
        shapes[k].tr.clear();
        shapes[k].rc.clear();
        for (int i = 0; i < tr_tmp.size(); ++i) {
            if (keep_tr[i]) {
                shapes[k].tr.push_back(tr_tmp[i]);
            }
        }
        for (int i = 0; i < rc_tmp.size(); ++i) {
            if (keep_rc[i]) {
                shapes[k].rc.push_back(rc_tmp[i]);
            }
        }
    }
    
    // output triangles
    cout << "static const face<3> triangles[" << triangles.size() << "] = {" << endl;
    for (int k = 0; k < triangles.size(); ++k) {
        cout << "    face<3> {" << endl;
        
        cout << "        { ";
        for (int i = 0; i < 3; ++i) {
            cout << (int)triangles[k].n.d[i] << ", ";
        }
        cout << (int)triangles[k].d << " }," << endl;
        cout << "        {" << endl;
        for (int i = 0; i < 3; ++i) {
            cout << "            { ";
            for (int j = 0; j < 3; ++j) {
                cout << (int)triangles[k].nn[i].d[j] << ", ";
            }
            cout << (int)triangles[k].dd[i] << " }";
            if (i < 2) {
                cout << ",";
            }
            cout << endl;
        }
        cout << "        }" << endl;
        cout << "    }";
        if (k < triangles.size()-1) {
            cout << ",";
        }
        cout << endl;
    }
    cout << "};" << endl << endl;
    
    // output rectangles
    cout << "static const face<4> rectangles["<< rectangles.size() << "] = {" << endl;
    for (int k = 0; k < rectangles.size(); ++k) {
        cout << "    face<4> {" << endl;
        
        cout << "        { ";
        for (int i = 0; i < 3; ++i) {
            cout << (int)rectangles[k].n.d[i] << ", ";
        }
        cout << (int)rectangles[k].d << " }," << endl;
        cout << "        {" << endl;
        for (int i = 0; i < 4; ++i) {
            cout << "            { ";
            for (int j = 0; j < 3; ++j) {
                cout << (int)rectangles[k].nn[i].d[j] << ", ";
            }
            cout << (int)rectangles[k].dd[i] << " }";
            if (i < 3) {
                cout << ",";
            }
            cout << endl;
        }
        cout << "        }" << endl;
        cout << "    }";
        if (k < rectangles.size()-1) {
            cout << ",";
        }
        cout << endl;
    }
    cout << "};" << endl << endl;
    
    // output shapes
    cout << "static const shape shapes[256] = {" << endl;
    for (int k = 0; k < 256; ++k) {
        cout << "    shape {" << endl;
        cout << "        {";
        for (int i = 0; i < shapes[k].tr.size(); ++i) {
            if (i > 0) {
                cout << ", ";
            } else {
                cout << " ";
            }
            cout << shapes[k].tr[i];
        }
        cout << " }," << endl;
        cout << "        {";
        for (int i = 0; i < shapes[k].rc.size(); ++i) {
            if (i > 0) {
                cout << ", ";
            } else {
                cout << " ";
            }
            cout << shapes[k].rc[i];
        }
        cout << " }" << endl;
        cout << "    }";
        if (k < 255) {
            cout << ",";
        }
        cout << endl;
    }
    cout << "};" << endl;
    return 0;
}
