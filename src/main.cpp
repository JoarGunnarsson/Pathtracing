#include <iostream>
#include <fstream>
#include <chrono>
#include <stdexcept>
#include <thread>
#include "vec3.h"
#include "colors.h"
#include "denoise.h"
#include "objects.h"
#include "camera.h"
#include "utils.h"
#include "constants.h"
#include "objectunion.h"
#include "scene.h"

std::mutex progress_mutex;
std::mutex job_mutex;

struct PixelData {
    vec3 pixel_color = vec3(0, 0, 0);
    vec3 pixel_position = vec3(0, 0, 0);
    vec3 pixel_normal = vec3(0, 0, 0);
};

struct Job {
    int start_idx;
    int number_of_pixels;
    int number_of_samples;
};

struct ThreadContext {
    int* thread_progress;
    int thread_idx;
    std::mutex* pixel_mutexes;
    int* pixel_counters;
};

PixelData raytrace(Ray ray, Object const* const* objects, const int number_of_objects, Medium* background_medium) {
    PixelData data;
    vec3 color = vec3(0, 0, 0);
    vec3 throughput = vec3(1, 1, 1);
    double random_threshold = 1;
    bool allow_recursion = true;
    bool has_hit_surface = false;

    vec3 saved_point;
    double scatter_pdf;

    Medium* medium = background_medium;

    for (int depth = 0; depth <= constants::max_recursion_depth; depth++) {
        if (!medium) {
            medium = background_medium;
        }
        double scatter_distance = medium->sample_distance();

        ray.t_max = scatter_distance;
        Hit ray_hit;
        if (!find_closest_hit(ray_hit, ray, objects, number_of_objects)) {
            if (scatter_distance == constants::max_ray_distance) {
                break;
            }
            ray_hit.distance = constants::max_ray_distance;
        }

        bool scatter = scatter_distance < ray_hit.distance;
        scatter_distance = std::min(scatter_distance, ray_hit.distance);
        if (scatter) {
            color += medium->sample_emission() * throughput;
        }

        throughput *= medium->sample(objects, number_of_objects, scatter_distance, scatter);

        if (scatter) {
            vec3 scatter_point = ray.starting_position + ray.direction_vector * scatter_distance;
            vec3 scattered_direction = medium->sample_direction(ray.direction_vector);
            if (constants::enable_next_event_estimation) {
                ray_hit.intersection_point = scatter_point;

                color +=
                    sample_light(ray_hit, objects, number_of_objects, background_medium, medium, true) * throughput;

                ray.type = DIFFUSE;
                scatter_pdf = medium->phase_function(ray.direction_vector, scattered_direction);
                saved_point = scatter_point;
            }

            ray.starting_position = scatter_point;
            ray.direction_vector = scattered_direction;
        }
        else {
            if (!has_hit_surface) {
                data.pixel_position = ray_hit.intersection_point;
                data.pixel_normal = ray_hit.normal_vector;
                has_hit_surface = true;
            }

            bool is_specular_ray = ray.type == REFLECTED || ray.type == TRANSMITTED;
            Object const* hit_object = objects[ray_hit.intersected_object_index];

            // If a light source is hit, compute the light_pdf based on the saved_point (previous hitpoint) and use MIS to add the light.
            // Could move this to a separate function, make it clearer what it is doing.
            if (hit_object->is_light_source()) {
                double weight;
                if (!constants::enable_next_event_estimation || depth == 0 || is_specular_ray) {
                    weight = 1;
                }
                else {
                    double light_pdf = objects[ray_hit.intersected_object_index]->light_pdf(
                        ray_hit.intersection_point, saved_point, ray_hit.primitive_ID);
                    weight = mis_weight(1, scatter_pdf, 1, light_pdf);
                }
                vec3 light_emittance = hit_object->get_light_emittance(ray_hit);

                color += weight * light_emittance * throughput;
            }

            if (constants::enable_next_event_estimation) {
                color +=
                    sample_light(ray_hit, objects, number_of_objects, background_medium, medium, false) * throughput;
            }

            BrdfData brdf_result = hit_object->sample(ray_hit);
            // TODO: Rename allow_direct_light!
            // TODO: Rename is_virtual_surface variable...
            bool is_virtual_surface =
                hit_object->get_material(ray_hit.primitive_ID)
                    ->allow_direct_light(); //This deviates from usual pattern of object method calling material method, but is better?
            if (is_virtual_surface) {
                brdf_result.type = ray.type;
            }
            else {
                scatter_pdf = brdf_result.pdf;
                saved_point = ray_hit.intersection_point;
            }
            throughput *= brdf_result.brdf_over_pdf;

            double incoming_dot_normal = dot_vectors(ray_hit.incident_vector, ray_hit.normal_vector);
            double outgoing_dot_normal = dot_vectors(brdf_result.outgoing_vector, ray_hit.normal_vector);

            bool penetrating_boundary = incoming_dot_normal * outgoing_dot_normal > 0;

            // TODO: Do the below part before sampling, so we can get the correct medium for refractive index etc?
            // TODO: Can save current_medium and next_medium, and pass that into sample and compute_direct_light.
            // TODO: Each material should have inside and outside media

            if (penetrating_boundary) {
                if (ray_hit.outside) {
                    medium = hit_object->get_material(ray_hit.primitive_ID)->internal_medium;
                }
                else {
                    medium = hit_object->get_material(ray_hit.primitive_ID)->external_medium;
                }
            }
            ray.starting_position = ray_hit.intersection_point;
            ray.direction_vector = brdf_result.outgoing_vector;
            ray.type = brdf_result.type;
        }

        if (depth < constants::min_recursion_steps) {
            random_threshold = 1;
            allow_recursion = true;
        }
        else {
            random_threshold = std::min(throughput.max(), 0.9);
            double random_value = random_uniform(0, 1);
            allow_recursion = random_value < random_threshold;
        }

        if (!allow_recursion) {
            break;
        }

        throughput /= random_threshold;
    }

    data.pixel_color = color;
    return data;
}

PixelData compute_pixel_color(const int x, const int y, const Scene& scene, const int number_of_samples) {
    PixelData data;
    vec3 pixel_color = vec3(0, 0, 0);
    for (int i = 0; i < number_of_samples; i++) {
        Ray ray;
        ray.starting_position = scene.camera.position;
        ray.type = TRANSMITTED;
        double new_x = x;
        double new_y = y;

        if (constants::enable_anti_aliasing) {
            new_x += random_normal() / 3.0;
            new_y += random_normal() / 3.0;
        }

        ray.direction_vector = scene.camera.get_starting_directions(new_x, new_y);
        PixelData sampled_data = raytrace(ray, scene.objects, scene.number_of_objects, scene.medium);
        data.pixel_position += sampled_data.pixel_position;
        data.pixel_normal += sampled_data.pixel_normal;
        pixel_color += sampled_data.pixel_color;
    }

    data.pixel_color = pixel_color / static_cast<double>(number_of_samples);
    data.pixel_position = data.pixel_position / static_cast<double>(number_of_samples);
    data.pixel_normal = data.pixel_normal / static_cast<double>(number_of_samples);
    return data;
}

void print_progress(double progress) {
    if (progress <= 1.0) {
        int bar_width = 60;

        std::clog << "[";
        int pos = static_cast<int>(bar_width * progress);
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos)
                std::clog << "=";
            else if (i == pos)
                std::clog << ">";
            else
                std::clog << " ";
        }
        std::clog << "] " << int(progress * 100.0) << " %\r";
    }
}

void update_and_report_progress(int* thread_progress) {
    int progress = 0;
    for (int i = 0; i < constants::number_of_threads; i++) {
        progress += thread_progress[i];
    }

    std::unique_lock<std::mutex> lock(progress_mutex);
    double total_pixel_count = static_cast<double>(constants::WIDTH * constants::HEIGHT * constants::samples_per_pixel);
    print_progress(static_cast<double>(progress) / static_cast<double>(total_pixel_count));
}

template<typename T>
T increment_average(const T previous, const T addition, const double samples_prev, const double samples_add) {
    return (previous * samples_prev + addition * samples_add) / (samples_prev + samples_add);
}

void raytrace_section(const Job& job, const Scene& scene, PixelBuffers& buffers, ThreadContext& thread_context) {
    int log_interval = job.number_of_pixels / 20;
    for (int i = 0; i < job.number_of_pixels; i++) {
        int idx = job.start_idx + i;

        int x = idx % constants::WIDTH;
        int y = constants::HEIGHT - idx / constants::WIDTH;

        PixelData data = compute_pixel_color(x, y, scene, job.number_of_samples);

        thread_context.pixel_mutexes[static_cast<size_t>(idx)].lock();
        double prev_samples = static_cast<double>(thread_context.pixel_counters[idx]);

        for (int j = 0; j < 3; j++) {
            buffers.image[3 * idx + j] =
                increment_average(buffers.image[3 * idx + j], data.pixel_color[j], prev_samples, job.number_of_samples);
        }
        buffers.position_buffer[idx] =
            increment_average(buffers.position_buffer[idx], data.pixel_position, prev_samples, job.number_of_samples);
        buffers.normal_buffer[idx] =
            increment_average(buffers.normal_buffer[idx], data.pixel_normal, prev_samples, job.number_of_samples);
        thread_context.pixel_counters[idx] += job.number_of_samples;

        thread_context.pixel_mutexes[static_cast<size_t>(idx)].unlock();

        if (i % log_interval == 0 && i != 0) {
            thread_context.thread_progress[thread_context.thread_idx] += log_interval * job.number_of_samples;
            update_and_report_progress(thread_context.thread_progress);
        }
    }
}

Job get_job(std::vector<Job>* job_queue) {
    std::unique_lock<std::mutex> lock(job_mutex);
    Job job = job_queue->back();
    job_queue->pop_back();
    return job;
}

void thread_job(std::vector<Job>* job_queue, const Scene& scene, PixelBuffers buffers, ThreadContext thread_context) {
    while (!job_queue->empty()) {
        Job job = get_job(job_queue);
        raytrace_section(job, scene, buffers, thread_context);
    }
}

void clear_scene(Scene& scene) {
    for (int i = 0; i < scene.number_of_objects; i++) {
        delete scene.objects[i];
    }

    delete[] scene.objects;
    delete scene.pointer_manager;
}

int main(int argc, char* argv[]) {
    std::chrono::steady_clock::time_point begin_build = std::chrono::steady_clock::now();
    if (argc != 3) {
        throw std::runtime_error(
            "Invalid arguments provided.\n"
            "Usage: main settings_file\n\n"

            "positional arguments:\n"
            "   scene_file                  scene file path, relative to main project directory.\n"
            "   settings_file               settings file path, relative to main project directory.\n");
    }
    load_settings(std::string(argv[2]));
    Scene scene = load_scene(std::string(argv[1]));

    std::chrono::steady_clock::time_point end_build = std::chrono::steady_clock::now();
    std::clog << "Time taken to build scene: "
              << std::chrono::duration_cast<std::chrono::seconds>(end_build - begin_build).count() << "[s]"
              << std::endl;

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    PixelBuffers pixel_buffers;
    int pixel_count = constants::WIDTH * constants::HEIGHT;
    size_t array_size = static_cast<size_t>(pixel_count);
    pixel_buffers.position_buffer = new vec3[array_size];
    pixel_buffers.normal_buffer = new vec3[array_size];

    size_t FILESIZE = static_cast<size_t>(pixel_count * 3) * sizeof(double);

    std::clog << "Running program with number of threads: " << constants::number_of_threads << ".\n";
    std::thread thread_array[constants::number_of_threads];

    int image_fd;
    pixel_buffers.image = create_mmap(constants::raw_file_name, FILESIZE, image_fd);

    print_progress(0);
    int pixels_per_thread = ceil_division(pixel_count, constants::number_of_threads);
    int* thread_progress = new int[static_cast<size_t>(constants::number_of_threads)];

    std::vector<Job> job_queue;

    int number_of_iterations = ceil_division(constants::samples_per_pixel, constants::samples_per_iteration);
    for (int j = 0; j < number_of_iterations; j++) {
        for (int i = 0; i < constants::number_of_threads; i++) {
            int start_idx = pixels_per_thread * i;
            int pixels_remaining = pixel_count - i * pixels_per_thread;
            int pixels_to_handle = std::min(pixels_per_thread, pixels_remaining);
            Job job;
            job.start_idx = start_idx;
            job.number_of_pixels = pixels_to_handle;
            job.number_of_samples = std::min(constants::samples_per_iteration,
                                             constants::samples_per_pixel - j * constants::samples_per_iteration);
            job_queue.push_back(job);
        }
    }

    std::mutex* pixel_mutexes = new std::mutex[static_cast<size_t>(pixel_count)];
    int* pixel_counters = new int[static_cast<size_t>(pixel_count)];
    for (int i = 0; i < constants::number_of_threads; i++) {
        ThreadContext thread_context{thread_progress, i, pixel_mutexes, pixel_counters};
        thread_array[i] = std::thread(thread_job, &job_queue, scene, pixel_buffers, thread_context);
    }

    for (int i = 0; i < constants::number_of_threads; i++) {
        thread_array[i].join();
    }

    print_progress(1);
    std::clog << std::endl;

    if (constants::enable_atrous_filtering || constants::enable_median_filtering) {
        PixelBuffers denoising_buffers;
        int denoised_image_fd;
        denoising_buffers.image = create_mmap(constants::raw_denoised_file_name, FILESIZE, denoised_image_fd);

        std::memcpy(denoising_buffers.image, pixel_buffers.image, FILESIZE);
        denoising_buffers.position_buffer = pixel_buffers.position_buffer;
        denoising_buffers.normal_buffer = pixel_buffers.normal_buffer;

        denoise(denoising_buffers);
        close_mmap(denoising_buffers.image, FILESIZE, denoised_image_fd);
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
    std::clog << "Time taken: " << std::chrono::duration_cast<std::chrono::seconds>(end - begin).count() << "[s]"
              << std::endl;

    close_mmap(pixel_buffers.image, FILESIZE, image_fd);

    clear_scene(scene);

    delete[] pixel_buffers.position_buffer;
    delete[] pixel_buffers.normal_buffer;
    return 0;
}
