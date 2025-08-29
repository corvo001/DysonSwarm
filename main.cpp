// Dyson Swarm mejorado:
// - Órbitas elípticas (e, argp) y distribución Walker
// - HUD con potencia total/media y controles
// - LOD dinámico del tamaño de punto
// - Satélites color verde viridiano
// Controles: Ratón (LMB rotar, rueda zoom), P/O pausa, R reseed, +/- cambia N, F11 fullscreen sin bordes, Esc salir

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include <cmath>

#include <imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

// ------------------- Shaders -------------------
static const char* VS = R"(
#version 330 core
layout(location=0) in vec3 aPos;
layout(location=1) in vec3 iColor;
layout(location=2) in vec3 iPos;
uniform mat4 uMVP;
uniform float uPointSize;
out vec3 vColor;
void main(){
    gl_Position = uMVP * vec4(iPos + aPos, 1.0);
    gl_PointSize = uPointSize;
    vColor = iColor;
}
)";

static const char* FS = R"(
#version 330 core
in vec3 vColor;
out vec4 FragColor;
void main(){
    vec2 pc = gl_PointCoord*2.0 - 1.0;
    float d = dot(pc, pc);
    if(d > 1.0) discard;
    float edge = smoothstep(1.0, 0.8, 1.0 - d);
    FragColor = vec4(vColor, edge);
}
)";

static const char* VS_STAR = R"(
#version 330 core
layout(location=0) in vec3 aPos;
uniform mat4 uMVP;
void main(){
    gl_Position = uMVP * vec4(aPos, 1.0);
    gl_PointSize = 20.0;
}
)";

static const char* FS_STAR = R"(
#version 330 core
out vec4 FragColor;
void main(){
    vec2 pc = gl_PointCoord*2.0 - 1.0;
    float r2 = dot(pc, pc);
    if(r2>1.0) discard;
    float glow = exp(-3.0*r2);
    FragColor = vec4(1.0, 0.95, 0.6, 1.0) * glow;
}
)";

// ------------------- Estructuras -------------------
struct Body {
    float a, e, incl, lan, argp, w, phase;
    glm::vec3 color;
};

struct Camera {
    float distance = 10.0f;
    float yaw = 0.6f, pitch = 0.6f;
    double lastX=0, lastY=0;
    bool rotating=false;
    glm::mat4 view() const {
        float x = distance * cosf(pitch) * cosf(yaw);
        float y = distance * sinf(pitch);
        float z = distance * cosf(pitch) * sinf(yaw);
        return glm::lookAt(glm::vec3(x,y,z), glm::vec3(0,0,0), glm::vec3(0,1,0));
    }
};

// ------------------- Utilidades GL -------------------
static GLuint compile(GLenum type, const char* src){
    GLuint s = glCreateShader(type);
    glShaderSource(s, 1, &src, nullptr);
    glCompileShader(s);
    GLint ok; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
    if(!ok){
        GLint len=0; glGetShaderiv(s, GL_INFO_LOG_LENGTH, &len);
        std::string log(len, '\0');
        glGetShaderInfoLog(s, len, nullptr, log.data());
        std::cerr << "Shader error:\n" << log << std::endl;
        std::exit(1);
    }
    return s;
}

static GLuint makeProgram(const char* vsrc, const char* fsrc){
    GLuint vs = compile(GL_VERTEX_SHADER, vsrc);
    GLuint fs = compile(GL_FRAGMENT_SHADER, fsrc);
    GLuint prog = glCreateProgram();
    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glLinkProgram(prog);
    GLint ok; glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if(!ok){
        GLint len=0; glGetProgramiv(prog, GL_INFO_LOG_LENGTH, &len);
        std::string log(len, '\0');
        glGetProgramInfoLog(prog, len, nullptr, log.data());
        std::cerr << "Link error:\n" << log << std::endl;
        std::exit(1);
    }
    glDeleteShader(vs);
    glDeleteShader(fs);
    return prog;
}

struct SwarmGPU {
    GLuint vao=0, vbo_points=0, vbo_color=0, vbo_pos=0;
    size_t count=0;
    void init(){
        glGenVertexArrays(1, &vao);
        glBindVertexArray(vao);

        float origin[3] = {0,0,0};
        glGenBuffers(1, &vbo_points);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_points);
        glBufferData(GL_ARRAY_BUFFER, sizeof(origin), origin, GL_STATIC_DRAW);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);

        glGenBuffers(1, &vbo_color);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_color);
        glBufferData(GL_ARRAY_BUFFER, 1, nullptr, GL_DYNAMIC_DRAW);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
        glVertexAttribDivisor(1, 1);

        glGenBuffers(1, &vbo_pos);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_pos);
        glBufferData(GL_ARRAY_BUFFER, 1, nullptr, GL_DYNAMIC_DRAW);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
        glVertexAttribDivisor(2, 1);

        glBindVertexArray(0);
    }
    void update(size_t n, const std::vector<glm::vec3>& colors, const std::vector<glm::vec3>& pos){
        count = n;
        glBindVertexArray(vao);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_color);
        glBufferData(GL_ARRAY_BUFFER, colors.size()*sizeof(glm::vec3), colors.data(), GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_pos);
        glBufferData(GL_ARRAY_BUFFER, pos.size()*sizeof(glm::vec3), pos.data(), GL_DYNAMIC_DRAW);
        glBindVertexArray(0);
    }
    void draw(){
        glBindVertexArray(vao);
        glDrawArraysInstanced(GL_POINTS, 0, 1, (GLsizei)count);
        glBindVertexArray(0);
    }
};

// ------------------- Dinámica orbital -------------------
static float solve_kepler(float M, float e){
    float E = M;
    for(int i=0;i<6;++i){
        float f  = E - e*std::sin(E) - M;
        float fp = 1.0f - e*std::cos(E);
        E -= f/fp;
    }
    return E;
}

static glm::vec3 orbital_position(const Body& b, float t){
    float M = b.w * t + b.phase;
    M = std::fmod(M, glm::two_pi<float>());

    float e = glm::clamp(b.e, 0.0f, 0.9f);
    float E = solve_kepler(M, e);

    float cosE = std::cos(E), sinE = std::sin(E);
    float x_orb = b.a*(cosE - e);
    float y_orb = b.a*std::sqrt(1 - e*e)*sinE;

    float cO = std::cos(b.lan), sO = std::sin(b.lan);
    float ci = std::cos(b.incl), si = std::sin(b.incl);
    float cw = std::cos(b.argp), sw = std::sin(b.argp);

    // Rz(argp)
    float x1 =  x_orb*cw - y_orb*sw;
    float y1 =  x_orb*sw + y_orb*cw;
    float z1 =  0.0f;
    // Rx(incl)
    float x2 = x1;
    float y2 = y1*ci - z1*si;
    float z2 = y1*si + z1*ci;
    // Rz(LAN)
    float x3 = x2*cO - y2*sO;
    float y3 = x2*sO + y2*cO;
    float z3 = z2;

    return glm::vec3(x3,y3,z3);
}

// ------------------- Distribución inicial (Walker) -------------------
static void seed_swarm(std::vector<Body>& bodies, int N, uint32_t seed=0){
    if(seed==0) seed = (uint32_t)std::chrono::high_resolution_clock::now().time_since_epoch().count();
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist_r(2.5f, 6.0f);
    std::uniform_real_distribution<float> dist_e(0.0f, 0.2f);
    std::uniform_real_distribution<float> dist_argp(0.0f, glm::two_pi<float>());

    const int Q = std::max(1, (int)std::round(std::cbrt((double)N)));
    const int P = std::max(1, N / Q);
    const float inc0 = glm::radians(55.0f);
    const glm::vec3 viridian(0.25f, 0.51f, 0.43f);

    bodies.clear(); bodies.reserve(N);
    int k=0;
    for(int q=0; q<Q; ++q){
        float lan = (glm::two_pi<float>() * q) / Q;
        for(int p=0; p<P && k<N; ++p, ++k){
            Body b;
            b.a = dist_r(rng);
            b.e = dist_e(rng);
            b.incl = inc0;
            b.lan = lan;
            b.argp = dist_argp(rng);
            b.w = std::sqrt(1.0f/(b.a*b.a*b.a));   // μ=1
            b.phase = glm::two_pi<float>() * (float)p / (float)P;
            b.color = viridian;
            bodies.push_back(b);
        }
    }
}

// ------------------- Potencia -------------------
static float collector_power(const glm::vec3& pos){
    float r2 = glm::dot(pos,pos);
    if(r2<1e-6f) return 0.0f;
    float invr2 = 1.0f/r2;
    float shade = pos.y>0 ? 1.0f : 0.6f;               // “sombra” simple
    float sat = 1.0f - std::exp(-3.0f*invr2);          // saturación
    return shade * sat;
}

// ------------------- Fullscreen borderless toggle -------------------
struct WindowState {
    bool fullscreen = false;
    int prevX=100, prevY=100, prevW=1280, prevH=720;
};

static void toggle_borderless_fullscreen(GLFWwindow* win, WindowState& ws){
    if(!ws.fullscreen){
        glfwGetWindowPos(win, &ws.prevX, &ws.prevY);
        glfwGetWindowSize(win, &ws.prevW, &ws.prevH);
        GLFWmonitor* mon = glfwGetPrimaryMonitor();
        const GLFWvidmode* vm = glfwGetVideoMode(mon);
        glfwSetWindowAttrib(win, GLFW_DECORATED, GLFW_FALSE);
        glfwSetWindowAttrib(win, GLFW_FLOATING, GLFW_TRUE);
        glfwSetWindowMonitor(win, nullptr, 0, 0, vm->width, vm->height, vm->refreshRate);
        ws.fullscreen = true;
    } else {
        glfwSetWindowAttrib(win, GLFW_DECORATED, GLFW_TRUE);
        glfwSetWindowAttrib(win, GLFW_FLOATING, GLFW_FALSE);
        glfwSetWindowMonitor(win, nullptr, ws.prevX, ws.prevY, ws.prevW, ws.prevH, 0);
        ws.fullscreen = false;
    }
}

// ------------------- Main -------------------
int main(){
    if(!glfwInit()){ std::cerr << "GLFW init failed\n"; return 1; }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR,3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR,3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    GLFWwindow* win = glfwCreateWindow(1280,720,"Dyson Swarm",nullptr,nullptr);
    if(!win){ std::cerr << "Window creation failed\n"; glfwTerminate(); return 1; }
    glfwMakeContextCurrent(win);
    glfwSwapInterval(1);

    if(glewInit()!=GLEW_OK){ std::cerr << "GLEW init failed\n"; return 1; }

    GLuint progSwarm = makeProgram(VS,FS);
    GLuint progStar  = makeProgram(VS_STAR,FS_STAR);
    GLint locMVP_Swarm = glGetUniformLocation(progSwarm,"uMVP");
    GLint locPointSize = glGetUniformLocation(progSwarm,"uPointSize");
    GLint locMVP_Star  = glGetUniformLocation(progStar,"uMVP");

    SwarmGPU swarm; swarm.init();

    GLuint vaoStar=0, vboStar=0;
    glGenVertexArrays(1,&vaoStar);
    glBindVertexArray(vaoStar);
    glGenBuffers(1,&vboStar);
    glBindBuffer(GL_ARRAY_BUFFER,vboStar);
    float starPos[3]={0,0,0};
    glBufferData(GL_ARRAY_BUFFER,sizeof(starPos),starPos,GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,0,(void*)0);
    glBindVertexArray(0);

    int N=5000;
    std::vector<Body> bodies; seed_swarm(bodies,N);
    std::vector<glm::vec3> colors(N), positions(N);
    for(int i=0;i<N;++i){ colors[i]=bodies[i].color; positions[i]=orbital_position(bodies[i],0.0f); }
    swarm.update(N,colors,positions);

    Camera cam;
    bool paused=false;
    WindowState ws;

    // Input cámara
    glfwSetWindowUserPointer(win, &cam);
    glfwSetCursorPosCallback(win, [](GLFWwindow* w, double x, double y){
        auto* c = (Camera*)glfwGetWindowUserPointer(w);
        if(c->rotating){
            float dx = float(x - c->lastX);
            float dy = float(y - c->lastY);
            c->yaw   -= dx * 0.005f;
            c->pitch += dy * 0.005f;
            c->pitch = glm::clamp(c->pitch, -1.54f, 1.54f);
        }
        c->lastX = x; c->lastY = y;
    });
    glfwSetMouseButtonCallback(win, [](GLFWwindow* w, int button, int action, int){
        auto* c = (Camera*)glfwGetWindowUserPointer(w);
        if(button==GLFW_MOUSE_BUTTON_LEFT){
            c->rotating = (action==GLFW_PRESS);
        }
    });
    glfwSetScrollCallback(win, [](GLFWwindow* w, double, double yoff){
        auto* c = (Camera*)glfwGetWindowUserPointer(w);
        c->distance *= (yoff>0 ? 0.9f : 1.1f);
        c->distance = glm::clamp(c->distance, 2.0f, 200.0f);
    });

    glEnable(GL_PROGRAM_POINT_SIZE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_DEPTH_TEST);

    auto t0 = std::chrono::high_resolution_clock::now();
    double hudAccum=0.0, hudLast=std::chrono::duration<double>(t0.time_since_epoch()).count();
    int hudFrames=0; double fps=0.0;

    // ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(win,true);
    ImGui_ImplOpenGL3_Init("#version 330");

    float pointSize=3.0f;
    bool showOverlay=true;
    bool f11Down=false;

    while(!glfwWindowShouldClose(win)){
        glfwPollEvents();

        // Teclado
        if(glfwGetKey(win,GLFW_KEY_P)==GLFW_PRESS) paused=true;
        if(glfwGetKey(win,GLFW_KEY_O)==GLFW_PRESS) paused=false;
        if(glfwGetKey(win,GLFW_KEY_R)==GLFW_PRESS){ seed_swarm(bodies,N); for(int i=0;i<N;++i) colors[i]=bodies[i].color; }
        bool f11Pressed = (glfwGetKey(win,GLFW_KEY_F11)==GLFW_PRESS);
        if(f11Pressed && !f11Down) toggle_borderless_fullscreen(win, ws);
        f11Down = f11Pressed;
        if(glfwGetKey(win,GLFW_KEY_KP_ADD)==GLFW_PRESS || glfwGetKey(win,GLFW_KEY_EQUAL)==GLFW_PRESS){
            int newN = std::min(N+500, 50000);
            if(newN!=N){
                N=newN; seed_swarm(bodies,N);
                colors.assign(N,{}); positions.assign(N,{});
                for(int i=0;i<N;++i){ colors[i]=bodies[i].color; positions[i]=orbital_position(bodies[i],0.0f); }
            }
        }
        if(glfwGetKey(win,GLFW_KEY_KP_SUBTRACT)==GLFW_PRESS || glfwGetKey(win,GLFW_KEY_MINUS)==GLFW_PRESS){
            int newN = std::max(N-500, 500);
            if(newN!=N){
                N=newN; seed_swarm(bodies,N);
                colors.assign(N,{}); positions.assign(N,{});
                for(int i=0;i<N;++i){ colors[i]=bodies[i].color; positions[i]=orbital_position(bodies[i],0.0f); }
            }
        }
        if(glfwGetKey(win,GLFW_KEY_ESCAPE)==GLFW_PRESS){ glfwSetWindowShouldClose(win,1); }

        // Tiempo
        auto t1 = std::chrono::high_resolution_clock::now();
        float dt = std::chrono::duration<float>(t1 - t0).count();

        // Update sim
        if(!paused){
            for(int i=0;i<N;++i){
                positions[i] = orbital_position(bodies[i], dt);
            }
        }
        swarm.update(N,colors,positions);

        // Potencia
        float powerSum=0.0f;
        for(int i=0;i<N;++i) powerSum += collector_power(positions[i]);
        float powerAvg = (N>0)? powerSum/N : 0.0f;

        // Render
        int w,h; glfwGetFramebufferSize(win,&w,&h);
        glViewport(0,0,w,h);
        glClearColor(0.005f,0.01f,0.02f,1.0f);
        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

        glm::mat4 P=glm::perspective(glm::radians(50.0f),(float)w/(float)h,0.01f,1000.0f);
        glm::mat4 V=cam.view();
        glm::mat4 M(1.0f);
        glm::mat4 MVP=P*V*M;

        // ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        if(showOverlay){
            ImGui::Begin("Dyson HUD",&showOverlay,ImGuiWindowFlags_AlwaysAutoResize);
            ImGui::Text("FPS: %.1f", fps);
            ImGui::Text("Satellites: %d", N);
            ImGui::Text("Power total: %.3f", powerSum);
            ImGui::Text("Power avg  : %.5f", powerAvg);
            ImGui::Checkbox("Pause", &paused);
            ImGui::SliderFloat("PointSize", &pointSize, 1.0f, 8.0f, "%.1f");
            if(ImGui::Button("Reseed (R)")){
                seed_swarm(bodies,N);
                for(int i=0;i<N;++i) colors[i]=bodies[i].color;
            }
            ImGui::End();
        }

        // Estrella
        glUseProgram(progStar);
        glUniformMatrix4fv(locMVP_Star,1,GL_FALSE,glm::value_ptr(MVP));
        glBindVertexArray(vaoStar);
        glDrawArrays(GL_POINTS,0,1);

        // LOD dinámico para satélites
        glm::vec3 eye = glm::vec3(glm::inverse(cam.view())[3]); // posición cámara
        double acc=0.0; int samples=0;
        for(int i=0;i<N; i+=std::max(1, N/512)){
            acc += glm::length(positions[i] - eye);
            ++samples;
        }
        float distAvg = (samples>0)? float(acc/samples) : cam.distance;
        float lodScale = glm::clamp(15.0f / (0.1f + distAvg), 0.5f, 4.0f);
        float pointSizeDynamic = pointSize * lodScale;

        // Enjambre
        glUseProgram(progSwarm);
        glUniformMatrix4fv(locMVP_Swarm,1,GL_FALSE,glm::value_ptr(MVP));
        glUniform1f(locPointSize, pointSizeDynamic);
        swarm.draw();

        // HUD en título (FPS)
        double now = std::chrono::duration<double>(t1.time_since_epoch()).count();
        hudFrames++; hudAccum += (now - hudLast); hudLast = now;
        if(hudAccum >= 1.0) {
            fps = hudFrames / hudAccum; hudFrames = 0; hudAccum = 0.0;
            char title[256];
            snprintf(title, sizeof(title),
                     "Dyson Swarm | FPS: %.1f | N: %d | %s",
                     fps, N, paused ? "Pausa" : "Run");
            glfwSetWindowTitle(win, title);
        }

        // ImGui render
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());

        glfwSwapBuffers(win);
    }

    // Shutdown
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glDeleteProgram(progSwarm);
    glDeleteProgram(progStar);
    glfwDestroyWindow(win);
    glfwTerminate();
    return 0;
}
