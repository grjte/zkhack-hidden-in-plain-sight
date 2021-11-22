#![allow(unused, unreachable_code, dead_code)]
use std::ops::Div;
use std::ops::{Mul, Sub};

use ark_bls12_381::{Fr, G1Affine};
use ark_ec::{AffineCurve, ProjectiveCurve};
use ark_ff::*;
use ark_poly::{
    univariate::DensePolynomial, EvaluationDomain, GeneralEvaluationDomain, Polynomial,
    UVPolynomial,
};
use ark_serialize::CanonicalDeserialize;
use hidden_in_plain_sight::{generate::kzg_commit, PUZZLE_DESCRIPTION};
use prompt::{puzzle, welcome};

fn read_cha_from_file() -> (Vec<G1Affine>, Vec<Vec<Fr>>, Fr, Fr, G1Affine, Fr, Fr) {
    use std::fs::File;
    use std::io::prelude::*;

    let mut file = File::open("challenge_data").unwrap();
    let mut bytes: Vec<u8> = vec![];
    file.read_to_end(&mut bytes).unwrap();

    let setup_bytes: Vec<u8> = bytes[0..98312].to_vec();
    let accts_bytes: Vec<u8> = bytes[98312..1130320].to_vec();
    let cha_1_bytes: Vec<u8> = bytes[1130320..1130352].to_vec();
    let cha_2_bytes: Vec<u8> = bytes[1130352..1130384].to_vec();
    let commt_bytes: Vec<u8> = bytes[1130384..1130480].to_vec();
    let opn_1_bytes: Vec<u8> = bytes[1130480..1130512].to_vec();
    let opn_2_bytes: Vec<u8> = bytes[1130512..1130544].to_vec();

    let setup = Vec::<G1Affine>::deserialize_unchecked(&setup_bytes[..]).unwrap();
    let accts = Vec::<Vec<Fr>>::deserialize_unchecked(&accts_bytes[..]).unwrap();
    let cha_1 = Fr::deserialize_unchecked(&cha_1_bytes[..]).unwrap();
    let cha_2 = Fr::deserialize_unchecked(&cha_2_bytes[..]).unwrap();
    let commt = G1Affine::deserialize_unchecked(&commt_bytes[..]).unwrap();
    let opn_1 = Fr::deserialize_unchecked(&opn_1_bytes[..]).unwrap();
    let opn_2 = Fr::deserialize_unchecked(&opn_2_bytes[..]).unwrap();

    (setup, accts, cha_1, cha_2, commt, opn_1, opn_2)
}

fn main() {
    welcome();
    puzzle(PUZZLE_DESCRIPTION);

    let (setup, accts, cha_1, cha_2, commt, opn_1, opn_2) = read_cha_from_file();

    // the # of accounts used to generate the puzzle
    let number_of_accts = 1000usize;
    // get the roots of unity (use the same domain used to generate the puzzle)
    let domain: GeneralEvaluationDomain<Fr> =
        GeneralEvaluationDomain::new(number_of_accts + 2).unwrap();

    // evaluate the vanishing polynomial at the two challenge points
    let van_1_eval = domain.evaluate_vanishing_polynomial(cha_1);
    let van_2_eval = domain.evaluate_vanishing_polynomial(cha_2);

    // Iterate through the accounts to find the solution polynomial derived from the account & get its commitment
    // let solution_commitment = run_acct_evals(
    //     setup.clone(),
    //     accts.clone(),
    //     cha_1,
    //     cha_2,
    //     commt,
    //     opn_1,
    //     opn_2,
    //     domain,
    //     van_1_eval,
    //     van_2_eval,
    // );
    // assert_eq!(solution_commitment, commt);

    // generate the blinded account solution from the found solution index in the accts
    println!("recipient address found at index {}", SOLUTION_INDEX);
    let solution_blinded_acct = get_blinded_acct_poly(
        domain,
        accts[SOLUTION_INDEX].clone(),
        cha_1,
        cha_2,
        opn_1,
        opn_2,
        van_1_eval,
        van_2_eval,
    );
    let solution_commt: G1Affine = kzg_commit(&solution_blinded_acct, &setup);
    assert_eq!(solution_commt, commt);
}

// this is the solution index found by the execution of run_acct_evals
const SOLUTION_INDEX: usize = 535;

// use the 2 challenges, 2 openings, evaluations of the vanishing polynomial at the 2 openings, and eval domain to
// get the blinded polynomial for the specified account acct
fn get_blinded_acct_poly(
    domain: GeneralEvaluationDomain<Fr>,
    acct: Vec<Fr>,
    cha_1: Fr,
    cha_2: Fr,
    opn_1: Fr,
    opn_2: Fr,
    van_1_eval: Fr,
    van_2_eval: Fr,
) -> DensePolynomial<Fr> {
    let acct_poly = DensePolynomial::from_coefficients_vec(domain.ifft(&acct));

    // use the challenges and openings to get the blinding factors
    let b_factors = get_blinding_factors(
        cha_1,
        cha_2,
        opn_1,
        opn_2,
        acct_poly.clone(),
        van_1_eval,
        van_2_eval,
    );

    // create the blinded account polynomial with the blinding factors
    let blinding_poly = DensePolynomial::from_coefficients_vec(b_factors);
    let blinded_acct_poly = acct_poly + blinding_poly.mul_by_vanishing_poly(domain);

    // return the blinded account poly
    (blinded_acct_poly)
}

// use the 2 challenges, 2 openings, and evaluations of the vanishing polynomial at the 2 openings to
// compute the 2 unknown blinding factors for this account polynomial
fn get_blinding_factors(
    cha_1: Fr,
    cha_2: Fr,
    opn_1: Fr,
    opn_2: Fr,
    acct_poly: DensePolynomial<Fr>,
    van_eval_1: Fr,
    van_eval_2: Fr,
) -> Vec<Fr> {
    // get openings for account polynomial
    let acct_eval_1 = acct_poly.evaluate(&cha_1);
    let acct_eval_2 = acct_poly.evaluate(&cha_2);

    // opening of blinding poly @ cha_1 = (opn_1 - acct_eval_1) / van_eval_1
    let num_1 = opn_1.sub(acct_eval_1);
    let opn_blinding_1 = num_1.div(van_eval_1);
    // opening of blinding poly @ cha_2 = (opn_2 - acct_eval_2) / van_eval_2
    let num_2 = opn_2.sub(acct_eval_2);
    let opn_blinding_2 = num_2.div(van_eval_2);

    // b_2 = (opn_blinding_1 - opn_blinding_2) / (cha_1 - cha_2)
    let opn_diff = opn_blinding_1.sub(opn_blinding_2);
    let b_2 = opn_diff.div(cha_1.sub(cha_2));
    // b_1 = opn_blinding_1 - b_2 * cha_1
    let b_1 = opn_blinding_1.sub(b_2.mul(cha_1));

    // return blinding factors
    (vec![b_1, b_2])
}

// evaluate every account at the challenge points and solve for the blinded polynomial until we find a matching commitment
// return the solution commitment
fn run_acct_evals(
    setup: Vec<G1Affine>,
    accts: Vec<Vec<Fr>>,
    cha_1: Fr,
    cha_2: Fr,
    commt: G1Affine,
    opn_1: Fr,
    opn_2: Fr,
    domain: GeneralEvaluationDomain<Fr>,
    van_1_eval: Fr,
    van_2_eval: Fr,
) -> G1Affine {
    // initialize the solution vector
    let mut solution_commitment = G1Affine::zero();

    // brute force it: iterate through the accounts
    // use the 2 challenges and 2 openings to find the 2 unknown blinding factors for each account
    // the create the commitment for that blinded polynomial of each account and see if it matches the proof commitment
    // if so, we have found our account
    for index in 0..accts.len() {
        println!("checking acct at index {:?}", index);
        let acct = &accts[index];

        let blinded_acct_poly = get_blinded_acct_poly(
            domain,
            acct.clone(),
            cha_1,
            cha_2,
            opn_1,
            opn_2,
            van_1_eval,
            van_2_eval,
        );
        let commitment: G1Affine = kzg_commit(&blinded_acct_poly, &setup);

        if commt == commitment {
            println!("solution found at index {}", index);
            solution_commitment = commitment;
            break;
        }
    }

    // return the solution commitment
    solution_commitment
}
