use rand::prelude::*;

fn main() {
    //defining initial values
    let size_of_the_square_matrix=10000;
    let coupling_between_neighboring_spins=0.44;
    let applied_field=0.02;

    //defining initial configuration
    let initial_configuration: Vec<Vec<i32>>=(0..size_of_the_square_matrix)
    .map(|_| (0..size_of_the_square_matrix).map(|_| if rand::random::<bool>() { 1 } else { -1 }).collect())
    .collect();

    let final_configuration=sweep(initial_configuration,size_of_the_square_matrix
        ,coupling_between_neighboring_spins,applied_field);
    
    println!("{}",final_configuration[0][0])   
}

//Sweep function that implements the Monte Carlo Method 
//employing the Metropolis-Hastings Algorithm
fn sweep(mut initial_configuration:Vec<Vec<i32>>, size:i32,
    coupling_between_neighboring_spins:f64, 
    applied_field:f64) -> Vec<Vec<i32>>{
        let indexes_x: Vec<i32>= Vec::from_iter(0..size-1);
        let indexes_y: Vec<i32>= Vec::from_iter(0..size-1);
        let random_number : f64 = rand::thread_rng().gen();
        for i in &indexes_x{
            for j in &indexes_y{
                let index_x: usize = (*i).try_into().unwrap();
                let index_y: usize = (*j).try_into().unwrap();
                let size2: usize=size.try_into().unwrap();

                let spin_up= initial_configuration[index_x]
                [if index_y==size2-1{0} else {index_y+1}];
                let spin_down= initial_configuration[index_x]
                [if index_y==0{size2-1} else {index_y-1}];
                let spin_left= initial_configuration
                [if index_x==size2-1{0} else {index_x+1}][index_y];
                let spin_right= initial_configuration
                [if index_x==0{size2-1} else {index_x+1}][index_y];

                let initial_spin=initial_configuration[index_x][index_y];
                let final_spin=initial_spin*(-1);

                let initial_energy=-coupling_between_neighboring_spins
                *((initial_spin*(spin_down+spin_left+spin_right+spin_up)) as f64)
                + applied_field*((spin_up+spin_down+spin_left+spin_right+initial_spin)
                 as f64);

                let final_energy=-coupling_between_neighboring_spins
                *((final_spin*(spin_down+spin_left+spin_right+spin_up)) as f64)
                + applied_field*((spin_up+spin_down+spin_left+spin_right+final_spin)
                as f64);

                let change_in_energy=initial_energy-final_energy;

                let probability_of_change=change_in_energy.exp();

                initial_configuration[index_x][index_y]
                =if probability_of_change>random_number{final_spin}
                else{initial_spin};
            }
            
        }
        initial_configuration

    }
