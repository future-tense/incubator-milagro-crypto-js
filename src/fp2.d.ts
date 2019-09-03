import {FpCommon} from "./fp-common";
import {FP} from "./fp";
import {BIG} from "./big";

export class FP2 extends FpCommon<FP2> {

    a: FP;
    b: FP;

    public constructor(c: FP2 | FP | BIG | number, d?: FP2 | FP | BIG | number);

    public getA(): BIG;

    public getB(): BIG;

    public set(c: FP, D: FP): void;

    public seta(c: FP): void;

    public bset(c: BIG, d: BIG): void;

    public bseta(c: BIG): void;

    public conj(): void;

    public div2(): void;

    public timesi(): void;

    public mul_ip(): void;

    public div_ip(): void;
}
