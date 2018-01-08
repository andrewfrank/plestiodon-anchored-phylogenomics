#!/bin/bash

# This list has loci with no information removed from it
LOCI=(L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13 L14 L15 L16 L17 L19 L21 L22 L23 L24 L26 L27 L28 L29 L30 L31 L32 L33 L34 L35 L36 L37 L38 L39 L40 L41 L42 L43 L44 L45 L46 L47 L48 L50 L51 L52 L53 L54 L55 L56 L57 L59 L60 L61 L62 L63 L65 L66 L67 L68 L69 L70 L71 L72 L73 L74 L75 L76 L77 L78 L79 L80 L81 L82 L83 L84 L86 L87 L88 L89 L90 L91 L92 L93 L94 L95 L96 L97 L98 L99 L100 L101 L102 L103 L105 L106 L107 L108 L109 L110 L111 L112 L113 L114 L115 L116 L118 L119 L120 L121 L122 L123 L124 L125 L126 L127 L128 L129 L130 L131 L132 L133 L134 L135 L136 L137 L138 L139 L140 L141 L142 L144 L145 L146 L147 L148 L149 L150 L151 L152 L153 L154 L155 L156 L158 L159 L160 L161 L162 L163 L164 L165 L166 L167 L168 L169 L173 L174 L175 L176 L177 L178 L179 L180 L181 L182 L183 L184 L185 L186 L187 L188 L189 L190 L191 L192 L193 L194 L195 L196 L197 L198 L199 L200 L201 L202 L205 L206 L207 L208 L209 L210 L211 L213 L214 L215 L216 L217 L218 L219 L220 L221 L222 L223 L224 L226 L227 L228 L229 L230 L231 L232 L233 L234 L235 L236 L237 L238 L239 L240 L241 L242 L243 L244 L245 L246 L247 L248 L249 L250 L251 L253 L254 L255 L256 L257 L258 L259 L260 L261 L262 L263 L264 L265 L266 L267 L268 L269 L270 L271 L272 L273 L274 L275 L276 L277 L278 L279 L280 L281 L282 L283 L284 L285 L286 L287 L288 L289 L290 L291 L292 L293 L294 L295 L296 L297 L298 L299 L300 L301 L302 L303 L304 L305 L306 L307 L308 L309 L310 L312 L313 L314 L315 L316 L317 L318 L319 L320 L321 L322 L323 L324 L325 L327 L328 L329 L330 L331 L332 L333 L334 L335 L336 L337 L338 L339 L341 L342 L343 L344 L345 L346 L347 L348 L349 L350 L351 L352 L353 L354 L355 L356 L357 L358 L359 L360 L361 L362 L363 L364 L365 L366 L367 L368 L369 L370 L371 L372 L373 L374 L375 L376 L377 L378 L379 L380 L381 L382 L383 L384 L385 L386 L387 L388 L389 L390 L391 L392 L393 L394 L395 L396 L397 L398 L399 L400 L401 L402 L403 L404 L405 L406 L407 L408 L409 L410 L411 L412)

INPUT_DIR="/Users/Andrew/Drive/Work/Projects/plestiodon-anchored-phylogenomics/data/original_seqs"
OUTPUT_DIR="/Users/Andrew/Drive/Work/Projects/plestiodon-anchored-phylogenomics/data/original_seqs_cleaned"

cd $INPUT_DIR
awk '/>/{sub(">","&"FILENAME"_");sub(/\T176_/,x);sub(/\.fasta/,x)}1' *.fasta > $OUTPUT_DIR/tmp.fasta

cd $OUTPUT_DIR
for LOCUS in "${LOCI[@]}"; do
  fgrep -v -e "I7494" -e "I18067" tmp.fasta |
  fgrep -A1 -e "$LOCUS"_ |
  sed '/^--$/d' |
  sed '/^$/d' > "$LOCUS".fasta
done

rm tmp.fasta
cat L*.fasta | paste - - | sort -k2,2 -g -tL | tr "\t" "\n" > concat-cleaned_seqs.fasta
