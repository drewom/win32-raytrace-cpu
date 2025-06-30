#include <windows.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <sal.h> // Include for annotations
#include <random>
#include <thread>

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

void RenderRaytraceToBitmap() {
    if (g_hBitmap) {
        DeleteObject(g_hBitmap);
        g_hBitmap = nullptr;
    }
    if (g_width <= 0 || g_height <= 0) return;

    float aspect = float(g_width) / g_height;
    float viewport_height = 2.0;
    float viewport_width = viewport_height * aspect;

    v3 lower_left(-viewport_width/2, -viewport_height/2, -1.0);
    v3 horizontal(viewport_width, 0.0, 0.0);
    v3 vertical(0.0, viewport_height, 0.0);
    v3 origin(0.0, 0.0, 0.0);

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
    auto render_band = [&](int y_start, int y_end) {
        for (int j = y_start; j < y_end; ++j) {
            for (int i = 0; i < g_width; ++i) {
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

                // Gamma correction (gamma=2.0)
                col = v3(std::sqrt(col.x), std::sqrt(col.y), std::sqrt(col.z));

                int ir = int(255.99f * (std::min)(1.0f, (std::max)(0.0f, col.x)));
                int ig = int(255.99f * (std::min)(1.0f, (std::max)(0.0f, col.y)));
                int ib = int(255.99f * (std::min)(1.0f, (std::max)(0.0f, col.z)));
                int idx = 4 * ((g_height-1-j)*g_width + i);
                pixels[idx+0] = (uint8_t)ib; // Blue
                pixels[idx+1] = (uint8_t)ig; // Green
                pixels[idx+2] = (uint8_t)ir; // Red
                pixels[idx+3] = 255;         // Alpha
            }
        }
    };

    int rows_per_thread = g_height / num_threads;
    int y = 0;
    for (unsigned int t = 0; t < num_threads; ++t) {
        int y_start = y;
        int y_end = (t == num_threads - 1) ? g_height : y + rows_per_thread;
        threads.emplace_back(render_band, y_start, y_end);
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