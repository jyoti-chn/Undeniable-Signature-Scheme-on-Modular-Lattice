# Undeniable-Signature-Scheme-on-Module-Lattice

## Overview
This project presents an undeniable signature scheme based on modular lattices, aiming to enhance security and align with proposed cryptographic standards. The GPV (Gentry, Peikert, and Vaikuntanathan) signature scheme serves as the base, with modifications to adapt it to modular lattices.

## Reasons for Modification
1. **Enhanced Security:** Shifting to modular lattices for improved security.
2. **Standardization Alignment:** Aligning with proposed cryptographic standards.
3. **Research Development:** Embracing algorithmic advances suited for modular lattices.

## Implementation Plan
Developed new algorithms for key generation, signing, and verification tailored for modular lattices.

## Key Algorithms

### KeyGen Algorithm
Generates signing and verification keys.
- Input: Security parameter n.
- Output: Signing key (sk) and verification key (pk).

### Sign Algorithm
Creates a signature for a given plaintext message.
- Input: Message (m) and signing key (sk).
- Output: Returns 1 for a valid signature; otherwise, 0.

### Confirmation Protocol
Interactive communication between Signer (S) and Verifier (V).
- Returns 1 for a valid signature; otherwise, 0.

### Disavowal Protocol
Interactive communication between Signer (S) and Verifier (V).
- Returns 1 for an invalid signature; otherwise, 0.

## Experimental Setup
Utilized two cryptographically secure collision-resistant hash functions:
- h: {0, 1}^ℤ → ℝ_𝑞^𝑙
- h1: {0, 1}𝑍 → ℝ_𝑞^(𝑙×𝑘)
Two variants for the modular hash function:
- Hash_com: To hash a random seed
- Hash_com_message: To hash a message concatenated with a random seed

## Key Generation
Executes MLTrapGen algorithm to obtain a pair (A, T).
- Outputs public key (PK) and secret key (SK).

## Signature Generation
Produces a signature σ on the message m.
- Outputs the signature σ = ("σ" _1, "σ" _2,"σ" _3) for the message m.

## Verification

2. **Confirmation/Disavowal Protocol:**
   This protocol operates between the signer and verifier 
   - If all verifications pass successfully, the verifier outputs 1; otherwise, it outputs 0.

