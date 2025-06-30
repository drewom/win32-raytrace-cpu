#include <windows.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <sal.h> // Include for annotations
#include <random>
#include <thread>
#include <emmintrin.h> // SIMD intrinsics
#include <immintrin.h> // For AVX512/SSE2 intrinsics
#include <chrono>

// Menu command IDs
#define IDM_AA_1X   1001
#define IDM_AA_2X   1002
#define IDM_AA_4X   1004
#define IDM_AA_8X   1008
#define IDM_AA_16X  1016

struct v3 {
    float x, y, z;
    v3(float x_=0, float y_=0, float z_=0) : x(x_), y(y_), z(z_) {}
    v3 operator+(const v3& b) const { return v3(x+b.x, y+b.y, z+b.z); }
    v3 operator-(const v3& b) const { return v3(x-b.x, y-b.y, z-b.z); }
    v3 operator*(float b) const { return v3(x*b, y*b, z*b); }
    v3 operator/(float b) const { return v3(x/b, y/b, z/b); }
    float dot(const v3& b) const { return x*b.x + y*b.y + z*b.z; }
    v3 norm() const { float mg=std::sqrtf(x*x+y*y+z*z); return *this/mg; }
};

struct ray {
    v3 orig, dir;
    ray(const v3& o, const v3& d) : orig(o), dir(d) {}
};

static bool hit_sphere(const v3& center, float radius, const ray& r, float& t) {
    v3 oc = r.orig - center;
    float a = r.dir.dot(r.dir);
    float b = 2.0f * oc.dot(r.dir);
    float c = oc.dot(oc) - radius*radius;
    float discriminant = b*b - 4*a*c;
    if (discriminant < 0) return false;
    t = (-b - std::sqrtf(discriminant)) / (2.0f*a);
    return t > 0;
}

static v3 ray_color(const ray& r) {
    float t;
    // Sphere at (0,0,-1) with radius 0.5
    if (hit_sphere(v3(0,0,-1), 0.5, r, t)) {
        v3 N = (r.orig + r.dir*t - v3(0,0,-1)).norm();
        return v3(N.x+1, N.y+1, N.z+1) * 0.5f;
    }
    // Simple ground plane at y = -0.5
    if (r.dir.y != 0) {
        t = (-0.5f - r.orig.y) / r.dir.y;
        if (t > 0) {
            v3 p = r.orig + r.dir*t;
            float checker = (int(std::floor(p.x) + std::floor(p.z)) % 2 == 0) ? 0.8f : 0.2f;
            return v3(checker, checker, checker);
        }
    }
    // Background gradient
    v3 unit_dir = r.dir.norm();
    t = 0.5f*(unit_dir.y + 1.0f);
    return v3(1.0f, 1.0f, 1.0f)*(1.0f-t) + v3(0.5f, 0.7f, 1.0f)*t;
}

// Globals for Win32
HBITMAP g_hBitmap = nullptr;
int g_width = 1000, g_height = 600;
int g_samples_per_pixel = 4; // Default AA

// Forward declaration
void RenderRaytraceToBitmap();

void UpdateMenuChecks(HMENU hMenu) {
    CheckMenuItem(hMenu, IDM_AA_1X,  MF_BYCOMMAND | (g_samples_per_pixel == 1  ? MF_CHECKED : MF_UNCHECKED));
    CheckMenuItem(hMenu, IDM_AA_2X,  MF_BYCOMMAND | (g_samples_per_pixel == 2  ? MF_CHECKED : MF_UNCHECKED));
    CheckMenuItem(hMenu, IDM_AA_4X,  MF_BYCOMMAND | (g_samples_per_pixel == 4  ? MF_CHECKED : MF_UNCHECKED));
    CheckMenuItem(hMenu, IDM_AA_8X,  MF_BYCOMMAND | (g_samples_per_pixel == 8  ? MF_CHECKED : MF_UNCHECKED));
    CheckMenuItem(hMenu, IDM_AA_16X, MF_BYCOMMAND | (g_samples_per_pixel == 16 ? MF_CHECKED : MF_UNCHECKED));
}

static LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
    switch (msg) {
    case WM_PAINT: {
        PAINTSTRUCT ps;
        HDC hdc = BeginPaint(hwnd, &ps);
        if (g_hBitmap) {
            HDC memDC = CreateCompatibleDC(hdc);
            HGDIOBJ oldbmp = SelectObject(memDC, g_hBitmap);
            BitBlt(hdc, 0, 0, g_width, g_height, memDC, 0, 0, SRCCOPY);
            SelectObject(memDC, oldbmp);
            DeleteDC(memDC);
        }
        EndPaint(hwnd, &ps);
        return 0;
    }
    case WM_SIZE: {
        int new_width = LOWORD(lParam);
        int new_height = HIWORD(lParam);
        if (new_width > 0 && new_height > 0) {
            g_width = new_width;
            g_height = new_height;
            RenderRaytraceToBitmap();
            InvalidateRect(hwnd, nullptr, FALSE);
        }
        return 0;
    }
    case WM_COMMAND: {
        HMENU hMenu = GetMenu(hwnd);
        int prev_samples = g_samples_per_pixel;
        switch (LOWORD(wParam)) {
        case IDM_AA_1X:   g_samples_per_pixel = 1; break;
        case IDM_AA_2X:   g_samples_per_pixel = 2; break;
        case IDM_AA_4X:   g_samples_per_pixel = 4; break;
        case IDM_AA_8X:   g_samples_per_pixel = 8; break;
        case IDM_AA_16X:  g_samples_per_pixel = 16; break;
        default: return DefWindowProc(hwnd, msg, wParam, lParam);
        }
        if (g_samples_per_pixel != prev_samples) {
            UpdateMenuChecks(hMenu);
            RenderRaytraceToBitmap();
            InvalidateRect(hwnd, nullptr, FALSE);
        }
        return 0;
    }
    case WM_DESTROY:
        if (g_hBitmap) DeleteObject(g_hBitmap);
        PostQuitMessage(0);
        return 0;
    default:
        return DefWindowProc(hwnd, msg, wParam, lParam);
    }
}

// Helper: Runtime check for AVX512F support
bool cpu_supports_avx512f() {
    int regs[4] = {0};
#if defined(_MSC_VER)
    __cpuidex(regs, 7, 0);
    return (regs[1] & (1 << 16)) != 0;
#elif defined(__GNUC__) || defined(__clang__)
    __asm__ __volatile__(
        "cpuid"
        : "=a"(regs[0]), "=b"(regs[1]), "=c"(regs[2]), "=d"(regs[3])
        : "a"(7), "c"(0));
    return (regs[1] & (1 << 16)) != 0;
#else
    return false;
#endif
}

void RenderRaytraceToBitmap() {
    using namespace std::chrono;
    auto start = high_resolution_clock::now();

    if (g_hBitmap) {
        DeleteObject(g_hBitmap);
        g_hBitmap = nullptr;
    }
    if (g_width <= 0 || g_height <= 0) return;

    float aspect = float(g_width) / g_height;
    float viewport_height = 2.0f;
    float viewport_width = viewport_height * aspect;

    v3 lower_left(-viewport_width/2, -viewport_height/2, -1.0f);
    v3 horizontal(viewport_width, 0.0f, 0.0f);
    v3 vertical(0.0f, viewport_height, 0.0f);
    v3 origin(0.0f, 0.0f, 0.0f);

    std::vector<uint8_t> pixels(g_width * g_height * 4);

    // Non-random subpixel supersampling (regular grid)
    int spp = g_samples_per_pixel;
    int grid_x = 1, grid_y = 1;
    for (int n = 1; n * n <= spp; ++n) {
        if (spp % n == 0) {
            grid_x = n;
            grid_y = spp / n;
        }
    }

    // Multithreading setup
    unsigned int num_threads = std::thread::hardware_concurrency();
    if (num_threads == 0) num_threads = 4; // Fallback
    std::vector<std::thread> threads;

    // Runtime AVX512F check (once per render)
    static bool use_avx512 = false;
#if defined(__AVX512F__)
        use_avx512 = cpu_supports_avx512f();
#endif 

    auto render_band_avx512 = [&](int y_start, int y_end) {
#if defined(__AVX512F__)
        for (int j = y_start; j < y_end; ++j) {
            int i = 0;
            for (; i <= g_width - 16; i += 16) {
                __m512 col_x = _mm512_setzero_ps();
                __m512 col_y = _mm512_setzero_ps();
                __m512 col_z = _mm512_setzero_ps();

                float r_x[16], r_y[16], r_z[16];

                for (int sy = 0; sy < grid_y; ++sy) {
                    for (int sx = 0; sx < grid_x; ++sx) {
                        float u[16], v[16];
                        for (int k = 0; k < 16; ++k) {
                            u[k] = float((i + k + (sx + 0.5f) / grid_x) / (g_width - 1));
                            v[k] = float((j + (sy + 0.5f) / grid_y) / (g_height - 1));
                        }
                        v3 dirs[16];
                        for (int k = 0; k < 16; ++k)
                            dirs[k] = lower_left + horizontal * u[k] + vertical * v[k];
                        for (int k = 0; k < 16; ++k) {
                            ray r(origin, dirs[k]);
                            v3 c = ray_color(r);
                            r_x[k] = c.x;
                            r_y[k] = c.y;
                            r_z[k] = c.z;
                        }
                        col_x = _mm512_add_ps(col_x, _mm512_loadu_ps(r_x));
                        col_y = _mm512_add_ps(col_y, _mm512_loadu_ps(r_y));
                        col_z = _mm512_add_ps(col_z, _mm512_loadu_ps(r_z));
                    }
                }
                float inv_spp = 1.0f / float(spp);
                __m512 inv = _mm512_set1_ps(inv_spp);
                col_x = _mm512_mul_ps(col_x, inv);
                col_y = _mm512_mul_ps(col_y, inv);
                col_z = _mm512_mul_ps(col_z, inv);

                // Gamma correction (sqrt)
                col_x = _mm512_sqrt_ps(col_x);
                col_y = _mm512_sqrt_ps(col_y);
                col_z = _mm512_sqrt_ps(col_z);

                // Clamp and convert to int
                __m512 max_val = _mm512_set1_ps(1.0f);
                __m512 min_val = _mm512_set1_ps(0.0f);
                col_x = _mm512_min_ps(_mm512_max_ps(col_x, min_val), max_val);
                col_y = _mm512_min_ps(_mm512_max_ps(col_y, min_val), max_val);
                col_z = _mm512_min_ps(_mm512_max_ps(col_z, min_val), max_val);

                __m512 scale = _mm512_set1_ps(255.99f);
                col_x = _mm512_mul_ps(col_x, scale);
                col_y = _mm512_mul_ps(col_y, scale);
                col_z = _mm512_mul_ps(col_z, scale);

                float ir[16], ig[16], ib[16];
                _mm512_storeu_ps(ir, col_x);
                _mm512_storeu_ps(ig, col_y);
                _mm512_storeu_ps(ib, col_z);

                for (int k = 0; k < 16; ++k) {
                    int idx = 4 * ((g_height - 1 - j) * g_width + (i + k));
                    pixels[idx + 0] = (uint8_t)ib[k]; // Blue
                    pixels[idx + 1] = (uint8_t)ig[k]; // Green
                    pixels[idx + 2] = (uint8_t)ir[k]; // Red
                    pixels[idx + 3] = 255;            // Alpha
                }
            }
            // Scalar fallback for remaining pixels
            for (; i < g_width; ++i) {
                v3 col(0, 0, 0);
                for (int sy = 0; sy < grid_y; ++sy) {
                    for (int sx = 0; sx < grid_x; ++sx) {
                        float u = (i + (sx + 0.5f) / grid_x) / (g_width-1);
                        float v = (j + (sy + 0.5f) / grid_y) / (g_height-1);
                        ray r(origin, lower_left + horizontal*u + vertical*v);
                        col = col + ray_color(r);
                    }
                }
                col = col / float(spp);
                col = v3(std::sqrtf(col.x), std::sqrtf(col.y), std::sqrtf(col.z));
                int ir = int(255.99f * (std::min)(1.0f, (std::max)(0.0f, col.x)));
                int ig = int(255.99f * (std::min)(1.0f, (std::max)(0.0f, col.y)));
                int ib = int(255.99f * (std::min)(1.0f, (std::max)(0.0f, col.z)));
                int idx = 4 * ((g_height-1-j)*g_width + i);
                pixels[idx+0] = (uint8_t)ib;
                pixels[idx+1] = (uint8_t)ig;
                pixels[idx+2] = (uint8_t)ir;
                pixels[idx+3] = 255;
            }
        }
#else
        // If not compiled with AVX512F, fallback to SSE2
        // (This branch will never be taken if __AVX512F__ is not defined)
#endif
    };

    auto render_band_sse2 = [&](int y_start, int y_end) {
        for (int j = y_start; j < y_end; ++j) {
            int i = 0;
            for (; i <= g_width - 4; i += 4) {
                __m128 col_x = _mm_setzero_ps();
                __m128 col_y = _mm_setzero_ps();
                __m128 col_z = _mm_setzero_ps();

                float r_x[4], r_y[4], r_z[4];

                for (int sy = 0; sy < grid_y; ++sy) {
                    for (int sx = 0; sx < grid_x; ++sx) {
                        float u[4], v[4];
                        for (int k = 0; k < 4; ++k) {
                            u[k] = float((i + k + (sx + 0.5f) / grid_x) / (g_width - 1));
                            v[k] = float((j + (sy + 0.5f) / grid_y) / (g_height - 1));
                        }
                        v3 dirs[4];
                        for (int k = 0; k < 4; ++k)
                            dirs[k] = lower_left + horizontal * u[k] + vertical * v[k];
                        for (int k = 0; k < 4; ++k) {
                            ray r(origin, dirs[k]);
                            v3 c = ray_color(r);
                            r_x[k] = c.x;
                            r_y[k] = c.y;
                            r_z[k] = c.z;
                        }
                        col_x = _mm_add_ps(col_x, _mm_loadu_ps(r_x));
                        col_y = _mm_add_ps(col_y, _mm_loadu_ps(r_y));
                        col_z = _mm_add_ps(col_z, _mm_loadu_ps(r_z));
                    }
                }
                float inv_spp = 1.0f / float(spp);
                __m128 inv = _mm_set1_ps(inv_spp);
                col_x = _mm_mul_ps(col_x, inv);
                col_y = _mm_mul_ps(col_y, inv);
                col_z = _mm_mul_ps(col_z, inv);

                // Gamma correction (sqrt)
                col_x = _mm_sqrt_ps(col_x);
                col_y = _mm_sqrt_ps(col_y);
                col_z = _mm_sqrt_ps(col_z);

                // Clamp and convert to int
                __m128 max_val = _mm_set1_ps(1.0f);
                __m128 min_val = _mm_set1_ps(0.0f);
                col_x = _mm_min_ps(_mm_max_ps(col_x, min_val), max_val);
                col_y = _mm_min_ps(_mm_max_ps(col_y, min_val), max_val);
                col_z = _mm_min_ps(_mm_max_ps(col_z, min_val), max_val);

                __m128 scale = _mm_set1_ps(255.99f);
                col_x = _mm_mul_ps(col_x, scale);
                col_y = _mm_mul_ps(col_y, scale);
                col_z = _mm_mul_ps(col_z, scale);

                float ir[4], ig[4], ib[4];
                _mm_storeu_ps(ir, col_x);
                _mm_storeu_ps(ig, col_y);
                _mm_storeu_ps(ib, col_z);

                for (int k = 0; k < 4; ++k) {
                    int idx = 4 * ((g_height - 1 - j) * g_width + (i + k));
                    pixels[idx + 0] = (uint8_t)ib[k]; // Blue
                    pixels[idx + 1] = (uint8_t)ig[k]; // Green
                    pixels[idx + 2] = (uint8_t)ir[k]; // Red
                    pixels[idx + 3] = 255;            // Alpha
                }
            }
            // Scalar fallback for remaining pixels
            for (; i < g_width; ++i) {
                v3 col(0, 0, 0);
                for (int sy = 0; sy < grid_y; ++sy) {
                    for (int sx = 0; sx < grid_x; ++sx) {
                        float u = (i + (sx + 0.5f) / grid_x) / (g_width-1);
                        float v = (j + (sy + 0.5f) / grid_y) / (g_height-1);
                        ray r(origin, lower_left + horizontal*u + vertical*v);
                        col = col + ray_color(r);
                    }
                }
                col = col / float(spp);
                col = v3(std::sqrtf(col.x), std::sqrtf(col.y), std::sqrtf(col.z));
                int ir = int(255.99f * (std::min)(1.0f, (std::max)(0.0f, col.x)));
                int ig = int(255.99f * (std::min)(1.0f, (std::max)(0.0f, col.y)));
                int ib = int(255.99f * (std::min)(1.0f, (std::max)(0.0f, col.z)));
                int idx = 4 * ((g_height-1-j)*g_width + i);
                pixels[idx+0] = (uint8_t)ib;
                pixels[idx+1] = (uint8_t)ig;
                pixels[idx+2] = (uint8_t)ir;
                pixels[idx+3] = 255;
            }
        }
    };

    // Launch threads
    int rows_per_thread = g_height / num_threads;
    int y = 0;
    for (unsigned int t = 0; t < num_threads; ++t) {
        int y_start = y;
        int y_end = (t == num_threads - 1) ? g_height : y + rows_per_thread;
        if (use_avx512) {
            threads.emplace_back(render_band_avx512, y_start, y_end);
        } else {
            threads.emplace_back(render_band_sse2, y_start, y_end);
        }
        y = y_end;
    }
    for (auto& th : threads) th.join();

    // Create DIB section
    BITMAPINFO bmi = {};
    bmi.bmiHeader.biSize = sizeof(BITMAPINFOHEADER);
    bmi.bmiHeader.biWidth = g_width;
    bmi.bmiHeader.biHeight = -g_height; // top-down
    bmi.bmiHeader.biPlanes = 1;
    bmi.bmiHeader.biBitCount = 32;
    bmi.bmiHeader.biCompression = BI_RGB;

    void* pBits = nullptr;
    g_hBitmap = CreateDIBSection(nullptr, &bmi, DIB_RGB_COLORS, &pBits, nullptr, 0);
    if (g_hBitmap && pBits) {
        memcpy(pBits, pixels.data(), pixels.size());
    }

    // Update the window title with render resolution, time, and architecture
    auto end = high_resolution_clock::now();
    auto ms = duration_cast<milliseconds>(end - start).count();
    HWND hwnd = FindWindow(L"RaytraceWindow", nullptr);
    if (hwnd) {
        wchar_t arch[32];
#if defined(__AVX512F__)
        if (use_avx512)
            wcscpy_s(arch, L"AVX512F");
        else
            wcscpy_s(arch, L"SSE2");
#else
        wcscpy_s(arch, L"SSE2");
#endif
        wchar_t title[160];
        swprintf(title, 160, L"Raytracing Demo - %dx%d - Last render: %lld ms [%s]", g_width, g_height, ms, arch);
        SetWindowText(hwnd, title);
    }
}

int WINAPI WinMain(_In_ HINSTANCE hInstance, _In_opt_ HINSTANCE hPrevInstance, _In_ LPSTR lpCmdLine, _In_ int nCmdShow) {
    // Register window class
    WNDCLASS wc = {};
    wc.lpfnWndProc = WndProc;
    wc.hInstance = hInstance;
    wc.lpszClassName = L"RaytraceWindow";
    wc.hCursor = LoadCursor(nullptr, IDC_ARROW);
    RegisterClass(&wc);

    // Create menu
    HMENU hMenu = CreateMenu();
    HMENU hAASub = CreateMenu();
    AppendMenu(hAASub, MF_STRING, IDM_AA_1X,  L"1x");
    AppendMenu(hAASub, MF_STRING, IDM_AA_2X,  L"2x");
    AppendMenu(hAASub, MF_STRING, IDM_AA_4X,  L"4x");
    AppendMenu(hAASub, MF_STRING, IDM_AA_8X,  L"8x");
    AppendMenu(hAASub, MF_STRING, IDM_AA_16X, L"16x");
    AppendMenu(hMenu, MF_POPUP, (UINT_PTR)hAASub, L"Antialiasing");

    HWND hwnd = CreateWindowEx(
        0, L"RaytraceWindow", L"Raytracing Demo",
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, g_width+16, g_height+60,
        nullptr, nullptr, hInstance, nullptr);

    SetMenu(hwnd, hMenu);
    UpdateMenuChecks(hMenu);

    ShowWindow(hwnd, nCmdShow);
    UpdateWindow(hwnd);

    // Initial render
    RenderRaytraceToBitmap();
    InvalidateRect(hwnd, nullptr, FALSE);

    // Message loop
    MSG msg;
    while (GetMessage(&msg, nullptr, 0, 0)) {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }
    return 0;
}