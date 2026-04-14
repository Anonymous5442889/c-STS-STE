use std::hint::black_box;

use c_ste::{dealer::Dealer, member::Member};
use criterion::{Criterion, criterion_group, criterion_main};

type E = ark_bls12_381::Bls12_381;

fn bench_hints(c: &mut Criterion) {
    let mut group = c.benchmark_group("single hint time");

    for k in 6..15 {
        let n = 1 << k;
        let (d, _) = Dealer::<E>::init_fast(n - 1);
        let members: Vec<Member<E>> = (0..n).map(|id| Member::<E>::new(id)).collect();
        
        if k <= 9 {
            group.bench_function(format!("Generate hint naive n = {}", n), |b| {
                b.iter(|| {
                    let hint = members[3].hint_gen_naive(&d);
                    let _ = black_box(hint);
                });
            });
        }
        
        group.bench_function(format!("Generate hint fast n = {}", n), |b| {
            b.iter(|| {
                let hint = members[3].hint_gen_lagrange(&d);
                let _ = black_box(hint);
            });
        });
    }

    group.finish();
}

criterion_group!(benches, bench_hints);
criterion_main!(benches);