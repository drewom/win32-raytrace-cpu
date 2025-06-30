#include <windows.h>
#include <cmath>
#include <vector>
#include <iostream>
#include <sal.h> // Include for annotations
#include <random>

struct v3 {
    double x, y, z;
    v3(double x_=0, double y_=0, double z_=0) : x(x_), y(y_), z(z_) {}
    v3 operator+(const v3& b) const { return v3(x+b.x, y+b.y, z+b.z); }
    v3 operator-(const v3& b) const { return v3(x-b.x, y-b.y, z-b.z); }
    v3 operator*(double b) const { return v3(x*b, y*b, z*b); }
    v3 operator/(double b) const { return v3(x/b, y/b, z/b); }
    double dot(const v3& b) const { return x*b.x + y*b.y + z*b.z; }
    v3 norm() const { double mg=std::sqrt(x*x+y*y+z*z); return *this/mg; }
};

struct ray {
    v3 orig, dir;
    ray(const v3& o, const v3& d) : orig(o), dir(d) {}
};

static bool hit_sphere(const v3& center, double radius, const ray& r, double& t) {
    v3 oc = r.orig - center;
    double a = r.dir.dot(r.dir);
    double b = 2.0 * oc.dot(r.dir);
    double c = oc.dot(oc) - radius*radius;
    double discriminant = b*b - 4*a*c;
    if (discriminant < 0) return false;
    t = (-b - std::sqrt(discriminant)) / (2.0*a);
    return t > 0;
}

static v3 ray_color(const ray& r) {
    double t;
    // Sphere at (0,0,-1) with radius 0.5
    if (hit_sphere(v3(0,0,-1), 0.5, r, t)) {
        v3 N = (r.orig + r.dir*t - v3(0,0,-1)).norm();
        return v3(N.x+1, N.y+1, N.z+1) * 0.5;
    }
    // Simple ground plane at y = -0.5
    if (r.dir.y != 0) {
        t = (-0.5 - r.orig.y) / r.dir.y;
        if (t > 0) {
            v3 p = r.orig + r.dir*t;
            double checker = (int(std::floor(p.x) + std::floor(p.z)) % 2 == 0) ? 0.8 : 0.2;
            return v3(checker, checker, checker);
        }
    }
    // Background gradient
    v3 unit_dir = r.dir.norm();
    t = 0.5*(unit_dir.y + 1.0);
    return v3(1.0, 1.0, 1.0)*(1.0-t) + v3(0.5, 0.7, 1.0)*t;
}

// Globals for Win32
HBITMAP g_hBitmap = nullptr;
int g_width = 1000, g_height = 600;

// Forward declaration
void RenderRaytraceToBitmap();

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

    // Aspect ratio correction
    double aspect = double(g_width) / g_height;
    double viewport_height = 2.0;
    double viewport_width = viewport_height * aspect;

    v3 lower_left(-viewport_width/2, -viewport_height/2, -1.0);
    v3 horizontal(viewport_width, 0.0, 0.0);
    v3 vertical(0.0, viewport_height, 0.0);
    v3 origin(0.0, 0.0, 0.0);

    std::vector<uint8_t> pixels(g_width * g_height * 4);

    // Antialiasing: number of samples per pixel
    const int samples_per_pixel = 4;
    std::mt19937 rng((unsigned int)time(nullptr));
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    for (int j = g_height-1; j >= 0; --j) {
        for (int i = 0; i < g_width; ++i) {
            v3 col(0, 0, 0);
            for (int s = 0; s < samples_per_pixel; ++s) {
                double u = (i + dist(rng)) / (g_width-1);
                double v = (j + dist(rng)) / (g_height-1);
                ray r(origin, lower_left + horizontal*u + vertical*v);
                col = col + ray_color(r);
            }
            col = col / double(samples_per_pixel);

            // Gamma correction (gamma=2.0)
            col = v3(std::sqrt(col.x), std::sqrt(col.y), std::sqrt(col.z));

            int ir = int(255.99 * min(1.0, max(0.0, col.x)));
            int ig = int(255.99 * min(1.0, max(0.0, col.y)));
            int ib = int(255.99 * min(1.0, max(0.0, col.z)));
            int idx = 4 * ((g_height-1-j)*g_width + i);
            pixels[idx+0] = (uint8_t)ib; // Blue
            pixels[idx+1] = (uint8_t)ig; // Green
            pixels[idx+2] = (uint8_t)ir; // Red
            pixels[idx+3] = 255;         // Alpha
        }
    }

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

    HWND hwnd = CreateWindowEx(
        0, L"RaytraceWindow", L"Raytracing Demo",
        WS_OVERLAPPEDWINDOW,
        CW_USEDEFAULT, CW_USEDEFAULT, g_width+16, g_height+39,
        nullptr, nullptr, hInstance, nullptr);

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