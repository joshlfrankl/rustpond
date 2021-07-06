// Rustpond -- A Rust port of Nanopond, a teeny tiny artificial life virtual machine
// Rust code copyright 2020 Joshua Franklin
// Original C implementation of Nanopond copyright 2005 Adam Ierymenko

// Rustpond is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

use clap::{App, Arg};
use rand::{Rng, SeedableRng};
use rand_core::RngCore;
use rand_xoshiro::Xoshiro256StarStar;
use std::fs::File;
use std::io::LineWriter;
use std::io::Write;

// Constants
const MUT_RATE: f64 = 0.00005;
const INFLOW_FREQUENCY: u64 = 100;
const INFLOW_SIZE: u64 = 4000;
const POND_SIZE_X: usize = 512;
const POND_SIZE_Y: usize = 512;
const POND_DEPTH: usize = 256;
const FAILED_KILL_PENALTY: u64 = 2;
const POND_DEPTH_BYTES: usize = POND_DEPTH / 2;
const REPORT_FREQUENCY: u64 = 5_000_000;
const PREFETCH_DEPTH: u64 = 8;
const POND_SIZE: usize = POND_SIZE_X * POND_SIZE_Y;

#[derive(Copy, Clone, Debug)]
enum Direction {
    UP,
    DOWN,
    LEFT,
    RIGHT,
}

#[derive(Clone)]
struct Cell {
    id: u64,
    parent_id: u64,
    lineage: u64,
    generation: u64,
    energy: u64,
    genome: [u8; POND_DEPTH_BYTES],
}

struct Pond {
    pond: Vec<Cell>,
    exec_buffer: [u8; POND_DEPTH_BYTES],
    stack: [usize; POND_DEPTH],
    rng: Xoshiro256StarStar, // Used for injection and randomness when executing a cell
    rng_secondary: Xoshiro256StarStar, // Used for deciding which cell to execute next; needed for reproducibility w/ PREFETCH_DEPTH
    cell_id_ctr: u64,
    clock: u64,
}

impl Pond {
    fn new(seed: u64) -> Pond {
        let pond = vec![
            Cell {
                id: 0,
                parent_id: 0,
                lineage: 0,
                generation: 0,
                energy: 0,
                genome: [0; POND_DEPTH_BYTES]
            };
            POND_SIZE_X * POND_SIZE_Y
        ];
        let exec_buffer = [0xFF; POND_DEPTH_BYTES];
        let rng: Xoshiro256StarStar;
        if seed == 0 {
            rng = Xoshiro256StarStar::seed_from_u64(rand::thread_rng().gen());
        } else {
            rng = Xoshiro256StarStar::seed_from_u64(seed);
        }

        let mut rng_secondary = rng.clone();
        rng_secondary.long_jump(); // Moves rng_secondary ahead of rng by 2^192 next_u64() calls; the generated sequences will never overlap

        let stack = [0usize; POND_DEPTH];
        let cell_id_ctr = 0u64;
        let clock = 0u64;
        return Pond {
            pond,
            exec_buffer,
            rng,
            rng_secondary,
            stack,
            cell_id_ctr,
            clock,
        };
    }

    fn calc_shift(ptr: usize) -> usize {
        (1 - (ptr % 2)) * 4
    }

    fn xy_to_idx(x: usize, y: usize) -> usize {
        return POND_SIZE_X * y + x;
    }
    
    fn reset_exec_buffer(&mut self) -> () {
        self.exec_buffer = [0xFF; POND_DEPTH_BYTES];
    }

    pub fn get_neighbor(idx: usize, dir: Direction) -> usize {
        let t = match dir {
            Direction::UP => (idx - POND_SIZE_X) % POND_SIZE,
            Direction::DOWN => (idx + POND_SIZE_X) % POND_SIZE,
            Direction::LEFT => (idx - 1) % POND_SIZE,
            Direction::RIGHT => (idx + 1) % POND_SIZE,
        };
        return t;
    }

    pub fn inject_cell(&mut self, cell_idx: usize) {
        self.pond[cell_idx].id = self.cell_id_ctr;
        self.pond[cell_idx].parent_id = 0;
        self.pond[cell_idx].lineage = self.cell_id_ctr;
        self.pond[cell_idx].generation = 0;
        self.pond[cell_idx].energy += INFLOW_SIZE;
        self.rng.fill_bytes(&mut self.pond[cell_idx].genome);
        self.cell_id_ctr += 1;
    }

    pub fn access_allowed(target: &Cell, guess: u8, sense: bool, rand: u8) -> bool {
        if target.energy == 0 || target.parent_id == 0 {
            return true;
        }
        let n_match = (!((target.genome[0] & 0xF0) ^ (guess & 0xF0))).count_ones() as u8;
        // True sense: more matches = more probable; used for positive interactions
        // False sense: more matches = less probable; used for negative interactions
        if sense {
            return (rand & 0x0F) <= n_match;
        } else {
            return (rand & 0x0F) >= n_match;
        }
    }

    pub fn run_sim(&mut self, run_length: u64, write_png: bool, out_path: &str) -> () {
        let mut exec_cell_buffer = [0usize; PREFETCH_DEPTH as usize];

        let mut img_ctr = 0u64; // Use an image counter so it's easy to merge in FFMPEG
        let mut idx: usize = self.rng_secondary.gen_range(0, POND_SIZE);
        let mut idx_in: usize = self.rng.gen_range(0, POND_SIZE);


        for i in 0..PREFETCH_DEPTH as usize {
            exec_cell_buffer[i] = idx;
            idx = self.rng_secondary.gen_range(0, POND_SIZE);
        }

        let tp = rayon::ThreadPoolBuilder::new()
            .num_threads(2)
            .build()
            .unwrap();
        
        let csv_file = File::create(format!("{}/output.csv", out_path)).unwrap();
        let mut csv_file = LineWriter::new(csv_file);
        csv_file.write_all(b"clock,total_energy,total_active_cells,total_viable_replicators,max_generations\n").unwrap();
                
        while self.clock != run_length {
            if self.clock % REPORT_FREQUENCY == 0 {
                let p = self.pond.clone();
                let c = self.clock;
                let csv_ref = &mut csv_file;
                tp.install(move || {
                    println!("Iteration: {}m", c / 1_000_000);
                    if write_png {
                        Pond::render_pond_png(&p, img_ctr, out_path);
                    }
                    Pond::print_info(&p, c, csv_ref);
                });
                if write_png {
                    img_ctr += 1;
                }
            }

            if self.clock % INFLOW_FREQUENCY == 0 {
                self.inject_cell(idx_in);
                idx_in = self.rng.gen_range(0, POND_SIZE);
            }

            let buff_idx = (self.clock % PREFETCH_DEPTH) as usize;
            let prefetch_idx = ((self.clock + (PREFETCH_DEPTH - 1)) % PREFETCH_DEPTH) as usize;

            unsafe { 
                // Manual prefetching since we know which cells are next,
                // but hardware prefetchers might not be able to figure it out.
                // Significantly improves perf on a Ryzen 3900X.
                let t_ptr: *const Cell = &self.pond[exec_cell_buffer[prefetch_idx]];
                let t_ptr = t_ptr as *const i8;
                std::arch::x86_64::_mm_prefetch(t_ptr, std::arch::x86_64::_MM_HINT_T0);
            }
            
            self.exec_cell(exec_cell_buffer[buff_idx]);
            idx = self.rng_secondary.gen_range(0, POND_SIZE);
            exec_cell_buffer[buff_idx] = idx;

            self.clock += 1;
        }
    }

    pub fn print_info(p: &Vec<Cell>, clock: u64, csv_handle: &mut LineWriter<File>) -> () {
        let mut total_energy = 0u64;
        let mut total_active_cells = 0u64;
        let mut total_viable_replicators = 0u64;
        let mut max_generation = 0u64;
        for c in p.iter() {
            if c.energy > 0 {
                total_energy += c.energy;
                total_active_cells += 1;
                if c.generation > 2 {
                    total_viable_replicators += 1;
                }
                if c.generation > max_generation {
                    max_generation = c.generation;
                }
            }
        }
        println!("Total energy: {}, Total active cells: {}, Total viable replicators: {}, Max generation: {}\n",
         total_energy, total_active_cells, total_viable_replicators, max_generation);
         csv_handle.write_all(format!("{},{},{},{},{}\n", clock, total_energy, total_active_cells,
                                    total_viable_replicators, max_generation).as_bytes()).unwrap();
    }


    pub fn exec_cell(&mut self, cell_idx: usize) -> () {
        if self.pond[cell_idx].energy == 0 {
            return;
        }
        let mut stop = false;
        let mut inst_ptr = 0usize;
        let mut mem_ptr = 0usize;
        let mut reg = 0u8;
        let mut facing = Direction::RIGHT;
        let mut loop_stack_ptr = 0;
        let mut false_loop_depth = 0;
        self.reset_exec_buffer();

        while !stop && self.pond[cell_idx].energy > 0 {
            self.pond[cell_idx].energy -= 1;
            let mut inst = (self.pond[cell_idx].genome[inst_ptr / 2] >> Pond::calc_shift(inst_ptr)) & 0x0F;

            if self.rng.gen_bool(MUT_RATE) {
                let rand_byte: u8 = self.rng.gen();
                if rand_byte & 0x80 > 0 {
                    inst = rand_byte & 0x0F;
                } else {
                    reg = rand_byte & 0x0F;
                }
            }

            if false_loop_depth > 0 {
                if inst == 0x09 {
                    false_loop_depth += 1;
                } else if inst == 0x0A {
                    false_loop_depth -= 1;
                }
            } else {
                match inst {
                    0x00 => {
                        // Reset VM status
                        reg = 0;
                        mem_ptr = 0;
                        facing = Direction::RIGHT;
                    }
                    0x01 => {
                        // Inc ptr
                        mem_ptr = (mem_ptr + 1) % POND_DEPTH;
                    } 
                    0x02 => {
                        // Dec ptr
                        mem_ptr = (mem_ptr - 1) % POND_DEPTH;
                    } 
                    0x03 => {
                        // Inc reg
                        reg = (reg + 1) & 0x0F;
                    } 
                    0x04 => {
                        // Dec Reg
                        reg = (reg - 1) & 0x0F;
                    } 
                    0x05 => {
                        // Read genome to reg
                        reg = self.pond[cell_idx].genome[mem_ptr / 2] >> Pond::calc_shift(mem_ptr) & 0x0F;
                    } 
                    0x06 => {
                        // Write reg to genome
                        self.pond[cell_idx].genome[mem_ptr / 2] &= !(0x0F << Pond::calc_shift(mem_ptr));
                        self.pond[cell_idx].genome[mem_ptr / 2] |= reg << Pond::calc_shift(mem_ptr);
                    }
                    0x07 => {
                        // Read buffer to reg
                        reg = self.exec_buffer[mem_ptr / 2] >> Pond::calc_shift(mem_ptr) & 0x0F;
                    } 
                    0x08 => {
                        // Write reg to buffer
                        self.exec_buffer[mem_ptr / 2] &= !(0x0F << Pond::calc_shift(mem_ptr));
                        self.exec_buffer[mem_ptr / 2] |= reg << Pond::calc_shift(mem_ptr);
                    }
                    0x09 => {
                        // Loop
                        if reg != 0 {
                            if loop_stack_ptr >= POND_DEPTH {
                                // Stack overflop; stop execution.
                                break;
                            } else {
                                self.stack[loop_stack_ptr] = inst_ptr;
                                loop_stack_ptr += 1;
                            }
                        } else {
                            false_loop_depth = 1;
                        }
                    }
                    0x0A => {
                        // Rep
                        if loop_stack_ptr > 0 {
                            loop_stack_ptr -= 1;
                            if reg != 0 {
                                inst_ptr = self.stack[loop_stack_ptr];
                                continue;
                            }
                        }
                    }
                    0x0B => {
                        // Turn
                        facing = match reg % 4 {
                            0 => Direction::UP,
                            1 => Direction::DOWN,
                            2 => Direction::LEFT,
                            3 => Direction::RIGHT,
                            _ => Direction::RIGHT,
                        };
                    }
                    0x0C => {
                        // Exchange
                        inst_ptr = (inst_ptr + 1) % POND_DEPTH;
                        let tmp = reg;
                        reg = self.pond[cell_idx].genome[inst_ptr / 2] >> Pond::calc_shift(inst_ptr) & 0x0F;
                        self.pond[cell_idx].genome[inst_ptr / 2] &= !(0x0F << Pond::calc_shift(inst_ptr));
                        self.pond[cell_idx].genome[inst_ptr / 2] |= tmp << Pond::calc_shift(inst_ptr);
                    }
                    0x0D => {
                        // Kill
                        let neighbor_idx = Pond::get_neighbor(cell_idx, facing);
                        if Pond::access_allowed(
                            &self.pond[neighbor_idx],
                            reg,
                            false,
                            self.rng.gen(),
                        ) {
                            self.pond[neighbor_idx].genome[0] = 0xFF;
                            self.pond[neighbor_idx].genome[1] = 0xFF;
                            self.pond[neighbor_idx].id = self.cell_id_ctr;
                            self.pond[neighbor_idx].parent_id = 0;
                            self.pond[neighbor_idx].lineage = self.cell_id_ctr;
                            self.pond[neighbor_idx].generation = 0;
                            self.cell_id_ctr += 1;
                        } else if self.pond[neighbor_idx].generation > 2 {
                            let penalty = self.pond[cell_idx].energy / FAILED_KILL_PENALTY;
                            let _ = self.pond[cell_idx].energy.saturating_sub(penalty);
                        }
                    }
                    0x0E => {
                        // Share; equalize energy between self and neighbor
                        let neighbor_idx = Pond::get_neighbor(cell_idx, facing);
                        if Pond::access_allowed(&self.pond[neighbor_idx], reg, true, self.rng.gen())
                        {
                            let total_energy =
                                self.pond[neighbor_idx].energy + self.pond[cell_idx].energy;
                            self.pond[cell_idx].energy = total_energy / 2;
                            self.pond[neighbor_idx].energy =
                                total_energy.saturating_sub(self.pond[cell_idx].energy);
                        }
                    }
                    0x0F => {
                        // Stop
                        stop = true;
                    } 
                    _ => return,
                }
            }
            inst_ptr = (inst_ptr + 1) % POND_DEPTH;
        }
        if self.exec_buffer[0] != 0xFF {
            let neighbor_idx = Pond::get_neighbor(cell_idx, facing);
            if (self.pond[neighbor_idx].energy > 0) && Pond::access_allowed(&self.pond[neighbor_idx], reg, false, self.rng.gen())
            {
                self.pond[neighbor_idx].id = self.cell_id_ctr;
                self.cell_id_ctr += 1;
                self.pond[neighbor_idx].parent_id = self.pond[cell_idx].id;
                self.pond[neighbor_idx].lineage = self.pond[cell_idx].lineage;
                self.pond[neighbor_idx].generation = self.pond[cell_idx].generation + 1;
                self.pond[neighbor_idx].genome = self.exec_buffer.clone();
            }
        }
    }

    pub fn render_pond_png(pond: &Vec<Cell>, id: u64, out_path: &str) -> () {
        let mut out_image = image::ImageBuffer::new(POND_SIZE_X as u32, POND_SIZE_Y as u32);
        let mut color;
        let mut c: &Cell;
        for (x, y, pixel) in out_image.enumerate_pixels_mut() {
            c = &pond[Pond::xy_to_idx(x as usize, y as usize)];
            if c.energy == 0 || c.generation < 2 {
                color = image::Rgb([0u8, 0u8, 0u8]);
            } else {
                let mut r: u32 = 0;
                for i in c.genome.iter() {
                    r += *i as u32;
                }
                color = image::Rgb([
                    (r & 0xF00 >> 2) as u8,
                    (r & 0x0F0 >> 1) as u8,
                    (r & 0x00F) as u8,
                ]);
            }
            *pixel = color;
        }
        out_image
            .save(format!("{}/rp_{}.png", out_path, id)) // 
            .unwrap(); // Will panic if unable to save the file for whatever reason (eg. no space, permissions, etc.)
    }

}

fn main() {
    let cmd_matches = App::new("Rustpond")
        .version("0.9.0")
        .author("Josh Franklin <joshlfrankl@gmail.com>")
        .about("A Rust port of Nanopond, a tiny artificial life VM. Original Nanopond C implementation written by Adam Ierymenko.")
        .arg(
            Arg::with_name("benchmark")
                .short("b")
                .long("benchmark")
                .help("Run a benchmark with a consistent seed for one billion updates; ignores all other flags"),
        )
        .arg(
            Arg::with_name("updates")
                .short("u")
                .long("updates")
                .takes_value(true)
                .help("Number of updates to run; default is 2^64 - 1"),
        )
        .arg(
            Arg::with_name("seed")
                .short("s")
                .long("seed")
                .takes_value(true)
                .help("Random seed to use; default is a randomly generated seed"),
        ).arg(
            Arg::with_name("out_path")
                .short("o")
                .long("out_path")
                .takes_value(true)
                .help("Output filepath to save PNGs; default is the current directory"),
        ).arg(
            Arg::with_name("write_png")
                .short("p")
                .long("write_pngs")
                .help("Writes PNG files every REPORT_FREQUENCY updates; default is false"),
        )
        .get_matches();

    if cmd_matches.is_present("benchmark") {
        println!("Running in benchmark mode");
        let mut pond = Pond::new(16481049608676553708); // A consistent seed for comparing across runs.
        let ti = std::time::Instant::now();
        pond.run_sim(1_000_000_000, false, ".");
        let elapsed = ti.elapsed();
        println!("Time taken for 1 billion cycles: {:?}", elapsed);
        return;
    }

    let mut num_updates = std::u64::MAX;
    if cmd_matches.is_present("updates") {
        num_updates = cmd_matches.value_of("updates").unwrap().parse().unwrap();
    }

    let mut seed = 0u64;
    if cmd_matches.is_present("seed") {
        seed = cmd_matches.value_of("seed").unwrap().parse().unwrap();
    }
    let mut pond = Pond::new(seed);
    let ti = std::time::Instant::now();
    
    let mut out_path = ".";
    if cmd_matches.is_present("out_path") {
        out_path = cmd_matches.value_of("out_path").unwrap();
    }

    let mut write_pngs = false;
    if cmd_matches.is_present("write_png") {
        write_pngs = true;
    }

    pond.run_sim(num_updates, write_pngs, out_path);
    let elapsed = ti.elapsed();
    println!("Time taken {:?}", elapsed);
}
